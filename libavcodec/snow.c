/*
 * Copyright (C) 2004 Michael Niedermayer <michaelni@gmx.at>
 *
 * This file is part of Libav.
 *
 * Libav is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * Libav is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Libav; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "libavutil/intmath.h"
#include "libavutil/log.h"
#include "libavutil/opt.h"
#include "avcodec.h"
#include "dsputil.h"
#include "dwt.h"
#include "snow.h"
#include "snowdata.h"

#include "rangecoder.h"
#include "mathops.h"
#include "h263.h"

#undef NDEBUG
#include <assert.h>


void ff_snow_inner_add_yblock(const uint8_t *obmc, const int obmc_stride, uint8_t * * block, int b_w, int b_h,
                              int src_x, int src_y, int src_stride, slice_buffer * sb, int add, uint8_t * dst8){
    int y, x;
    IDWTELEM * dst;
    for(y=0; y<b_h; y++){
        //FIXME ugly misuse of obmc_stride
        const uint8_t *obmc1= obmc + y*obmc_stride;
        const uint8_t *obmc2= obmc1+ (obmc_stride>>1);
        const uint8_t *obmc3= obmc1+ obmc_stride*(obmc_stride>>1);
        const uint8_t *obmc4= obmc3+ (obmc_stride>>1);
        dst = slice_buffer_get_line(sb, src_y + y);
        for(x=0; x<b_w; x++){
            int v=   obmc1[x] * block[3][x + y*src_stride]
                    +obmc2[x] * block[2][x + y*src_stride]
                    +obmc3[x] * block[1][x + y*src_stride]
                    +obmc4[x] * block[0][x + y*src_stride];

            v <<= 8 - LOG2_OBMC_MAX;
            if(FRAC_BITS != 8){
                v >>= 8 - FRAC_BITS;
            }
            if(add){
                v += dst[x + src_x];
                v = (v + (1<<(FRAC_BITS-1))) >> FRAC_BITS;
                if(v&(~255)) v= ~(v>>31);
                dst8[x + y*src_stride] = v;
            }else{
                dst[x + src_x] -= v;
            }
        }
    }
}

void ff_snow_reset_contexts(SnowContext *s){ //FIXME better initial contexts
    int plane_index, level, orientation;

    for(plane_index=0; plane_index<3; plane_index++){
        for(level=0; level<MAX_DECOMPOSITIONS; level++){
            for(orientation=level ? 1:0; orientation<4; orientation++){
                memset(s->plane[plane_index].band[level][orientation].state, MID_STATE, sizeof(s->plane[plane_index].band[level][orientation].state));
            }
        }
    }
    memset(s->header_state, MID_STATE, sizeof(s->header_state));
    memset(s->block_state, MID_STATE, sizeof(s->block_state));
}

int ff_snow_alloc_blocks(SnowContext *s){
    int w= -((-s->avctx->width )>>LOG2_MB_SIZE);
    int h= -((-s->avctx->height)>>LOG2_MB_SIZE);

    s->b_width = w;
    s->b_height= h;

    av_free(s->block);
    s->block= av_mallocz(w * h * sizeof(BlockNode) << (s->block_max_depth*2));
    return 0;
}

static void init_qexp(void){
    int i;
    double v=128;

    for(i=0; i<QROOT; i++){
        qexp[i]= lrintf(v);
        v *= pow(2, 1.0 / QROOT);
    }
}
static void mc_block(Plane *p, uint8_t *dst, const uint8_t *src, int stride, int b_w, int b_h, int dx, int dy){
    static const uint8_t weight[64]={
    8,7,6,5,4,3,2,1,
    7,7,0,0,0,0,0,1,
    6,0,6,0,0,0,2,0,
    5,0,0,5,0,3,0,0,
    4,0,0,0,4,0,0,0,
    3,0,0,5,0,3,0,0,
    2,0,6,0,0,0,2,0,
    1,7,0,0,0,0,0,1,
    };

    static const uint8_t brane[256]={
    0x00,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x11,0x12,0x12,0x12,0x12,0x12,0x12,0x12,
    0x04,0x05,0xcc,0xcc,0xcc,0xcc,0xcc,0x41,0x15,0x16,0xcc,0xcc,0xcc,0xcc,0xcc,0x52,
    0x04,0xcc,0x05,0xcc,0xcc,0xcc,0x41,0xcc,0x15,0xcc,0x16,0xcc,0xcc,0xcc,0x52,0xcc,
    0x04,0xcc,0xcc,0x05,0xcc,0x41,0xcc,0xcc,0x15,0xcc,0xcc,0x16,0xcc,0x52,0xcc,0xcc,
    0x04,0xcc,0xcc,0xcc,0x41,0xcc,0xcc,0xcc,0x15,0xcc,0xcc,0xcc,0x16,0xcc,0xcc,0xcc,
    0x04,0xcc,0xcc,0x41,0xcc,0x05,0xcc,0xcc,0x15,0xcc,0xcc,0x52,0xcc,0x16,0xcc,0xcc,
    0x04,0xcc,0x41,0xcc,0xcc,0xcc,0x05,0xcc,0x15,0xcc,0x52,0xcc,0xcc,0xcc,0x16,0xcc,
    0x04,0x41,0xcc,0xcc,0xcc,0xcc,0xcc,0x05,0x15,0x52,0xcc,0xcc,0xcc,0xcc,0xcc,0x16,
    0x44,0x45,0x45,0x45,0x45,0x45,0x45,0x45,0x55,0x56,0x56,0x56,0x56,0x56,0x56,0x56,
    0x48,0x49,0xcc,0xcc,0xcc,0xcc,0xcc,0x85,0x59,0x5A,0xcc,0xcc,0xcc,0xcc,0xcc,0x96,
    0x48,0xcc,0x49,0xcc,0xcc,0xcc,0x85,0xcc,0x59,0xcc,0x5A,0xcc,0xcc,0xcc,0x96,0xcc,
    0x48,0xcc,0xcc,0x49,0xcc,0x85,0xcc,0xcc,0x59,0xcc,0xcc,0x5A,0xcc,0x96,0xcc,0xcc,
    0x48,0xcc,0xcc,0xcc,0x49,0xcc,0xcc,0xcc,0x59,0xcc,0xcc,0xcc,0x96,0xcc,0xcc,0xcc,
    0x48,0xcc,0xcc,0x85,0xcc,0x49,0xcc,0xcc,0x59,0xcc,0xcc,0x96,0xcc,0x5A,0xcc,0xcc,
    0x48,0xcc,0x85,0xcc,0xcc,0xcc,0x49,0xcc,0x59,0xcc,0x96,0xcc,0xcc,0xcc,0x5A,0xcc,
    0x48,0x85,0xcc,0xcc,0xcc,0xcc,0xcc,0x49,0x59,0x96,0xcc,0xcc,0xcc,0xcc,0xcc,0x5A,
    };

    static const uint8_t needs[16]={
    0,1,0,0,
    2,4,2,0,
    0,1,0,0,
    15
    };

    int x, y, b, r, l;
    int16_t tmpIt   [64*(32+HTAPS_MAX)];
    uint8_t tmp2t[3][stride*(32+HTAPS_MAX)];
    int16_t *tmpI= tmpIt;
    uint8_t *tmp2= tmp2t[0];
    const uint8_t *hpel[11];
    assert(dx<16 && dy<16);
    r= brane[dx + 16*dy]&15;
    l= brane[dx + 16*dy]>>4;

    b= needs[l] | needs[r];
    if(p && !p->diag_mc)
        b= 15;

    if(b&5){
        for(y=0; y < b_h+HTAPS_MAX-1; y++){
            for(x=0; x < b_w; x++){
                int a_1=src[x + HTAPS_MAX/2-4];
                int a0= src[x + HTAPS_MAX/2-3];
                int a1= src[x + HTAPS_MAX/2-2];
                int a2= src[x + HTAPS_MAX/2-1];
                int a3= src[x + HTAPS_MAX/2+0];
                int a4= src[x + HTAPS_MAX/2+1];
                int a5= src[x + HTAPS_MAX/2+2];
                int a6= src[x + HTAPS_MAX/2+3];
                int am=0;
                if(!p || p->fast_mc){
                    am= 20*(a2+a3) - 5*(a1+a4) + (a0+a5);
                    tmpI[x]= am;
                    am= (am+16)>>5;
                }else{
                    am= p->hcoeff[0]*(a2+a3) + p->hcoeff[1]*(a1+a4) + p->hcoeff[2]*(a0+a5) + p->hcoeff[3]*(a_1+a6);
                    tmpI[x]= am;
                    am= (am+32)>>6;
                }

                if(am&(~255)) am= ~(am>>31);
                tmp2[x]= am;
            }
            tmpI+= 64;
            tmp2+= stride;
            src += stride;
        }
        src -= stride*y;
    }
    src += HTAPS_MAX/2 - 1;
    tmp2= tmp2t[1];

    if(b&2){
        for(y=0; y < b_h; y++){
            for(x=0; x < b_w+1; x++){
                int a_1=src[x + (HTAPS_MAX/2-4)*stride];
                int a0= src[x + (HTAPS_MAX/2-3)*stride];
                int a1= src[x + (HTAPS_MAX/2-2)*stride];
                int a2= src[x + (HTAPS_MAX/2-1)*stride];
                int a3= src[x + (HTAPS_MAX/2+0)*stride];
                int a4= src[x + (HTAPS_MAX/2+1)*stride];
                int a5= src[x + (HTAPS_MAX/2+2)*stride];
                int a6= src[x + (HTAPS_MAX/2+3)*stride];
                int am=0;
                if(!p || p->fast_mc)
                    am= (20*(a2+a3) - 5*(a1+a4) + (a0+a5) + 16)>>5;
                else
                    am= (p->hcoeff[0]*(a2+a3) + p->hcoeff[1]*(a1+a4) + p->hcoeff[2]*(a0+a5) + p->hcoeff[3]*(a_1+a6) + 32)>>6;

                if(am&(~255)) am= ~(am>>31);
                tmp2[x]= am;
            }
            src += stride;
            tmp2+= stride;
        }
        src -= stride*y;
    }
    src += stride*(HTAPS_MAX/2 - 1);
    tmp2= tmp2t[2];
    tmpI= tmpIt;
    if(b&4){
        for(y=0; y < b_h; y++){
            for(x=0; x < b_w; x++){
                int a_1=tmpI[x + (HTAPS_MAX/2-4)*64];
                int a0= tmpI[x + (HTAPS_MAX/2-3)*64];
                int a1= tmpI[x + (HTAPS_MAX/2-2)*64];
                int a2= tmpI[x + (HTAPS_MAX/2-1)*64];
                int a3= tmpI[x + (HTAPS_MAX/2+0)*64];
                int a4= tmpI[x + (HTAPS_MAX/2+1)*64];
                int a5= tmpI[x + (HTAPS_MAX/2+2)*64];
                int a6= tmpI[x + (HTAPS_MAX/2+3)*64];
                int am=0;
                if(!p || p->fast_mc)
                    am= (20*(a2+a3) - 5*(a1+a4) + (a0+a5) + 512)>>10;
                else
                    am= (p->hcoeff[0]*(a2+a3) + p->hcoeff[1]*(a1+a4) + p->hcoeff[2]*(a0+a5) + p->hcoeff[3]*(a_1+a6) + 2048)>>12;
                if(am&(~255)) am= ~(am>>31);
                tmp2[x]= am;
            }
            tmpI+= 64;
            tmp2+= stride;
        }
    }

    hpel[ 0]= src;
    hpel[ 1]= tmp2t[0] + stride*(HTAPS_MAX/2-1);
    hpel[ 2]= src + 1;

    hpel[ 4]= tmp2t[1];
    hpel[ 5]= tmp2t[2];
    hpel[ 6]= tmp2t[1] + 1;

    hpel[ 8]= src + stride;
    hpel[ 9]= hpel[1] + stride;
    hpel[10]= hpel[8] + 1;

    if(b==15){
        const uint8_t *src1= hpel[dx/8 + dy/8*4  ];
        const uint8_t *src2= hpel[dx/8 + dy/8*4+1];
        const uint8_t *src3= hpel[dx/8 + dy/8*4+4];
        const uint8_t *src4= hpel[dx/8 + dy/8*4+5];
        dx&=7;
        dy&=7;
        for(y=0; y < b_h; y++){
            for(x=0; x < b_w; x++){
                dst[x]= ((8-dx)*(8-dy)*src1[x] + dx*(8-dy)*src2[x]+
                         (8-dx)*   dy *src3[x] + dx*   dy *src4[x]+32)>>6;
            }
            src1+=stride;
            src2+=stride;
            src3+=stride;
            src4+=stride;
            dst +=stride;
        }
    }else{
        const uint8_t *src1= hpel[l];
        const uint8_t *src2= hpel[r];
        int a= weight[((dx&7) + (8*(dy&7)))];
        int b= 8-a;
        for(y=0; y < b_h; y++){
            for(x=0; x < b_w; x++){
                dst[x]= (a*src1[x] + b*src2[x] + 4)>>3;
            }
            src1+=stride;
            src2+=stride;
            dst +=stride;
        }
    }
}

void ff_snow_pred_block(SnowContext *s, uint8_t *dst, uint8_t *tmp, int stride, int sx, int sy, int b_w, int b_h, BlockNode *block, int plane_index, int w, int h){
    if(block->type & BLOCK_INTRA){
        int x, y;
        const unsigned color  = block->color[plane_index];
        const unsigned color4 = color*0x01010101;
        if(b_w==32){
            for(y=0; y < b_h; y++){
                *(uint32_t*)&dst[0 + y*stride]= color4;
                *(uint32_t*)&dst[4 + y*stride]= color4;
                *(uint32_t*)&dst[8 + y*stride]= color4;
                *(uint32_t*)&dst[12+ y*stride]= color4;
                *(uint32_t*)&dst[16+ y*stride]= color4;
                *(uint32_t*)&dst[20+ y*stride]= color4;
                *(uint32_t*)&dst[24+ y*stride]= color4;
                *(uint32_t*)&dst[28+ y*stride]= color4;
            }
        }else if(b_w==16){
            for(y=0; y < b_h; y++){
                *(uint32_t*)&dst[0 + y*stride]= color4;
                *(uint32_t*)&dst[4 + y*stride]= color4;
                *(uint32_t*)&dst[8 + y*stride]= color4;
                *(uint32_t*)&dst[12+ y*stride]= color4;
            }
        }else if(b_w==8){
            for(y=0; y < b_h; y++){
                *(uint32_t*)&dst[0 + y*stride]= color4;
                *(uint32_t*)&dst[4 + y*stride]= color4;
            }
        }else if(b_w==4){
            for(y=0; y < b_h; y++){
                *(uint32_t*)&dst[0 + y*stride]= color4;
            }
        }else{
            for(y=0; y < b_h; y++){
                for(x=0; x < b_w; x++){
                    dst[x + y*stride]= color;
                }
            }
        }
    }else{
        uint8_t *src= s->last_picture[block->ref].data[plane_index];
        const int scale= plane_index ?  s->mv_scale : 2*s->mv_scale;
        int mx= block->mx*scale;
        int my= block->my*scale;
        const int dx= mx&15;
        const int dy= my&15;
        const int tab_index= 3 - (b_w>>2) + (b_w>>4);
        sx += (mx>>4) - (HTAPS_MAX/2-1);
        sy += (my>>4) - (HTAPS_MAX/2-1);
        src += sx + sy*stride;
        if(   (unsigned)sx >= w - b_w - (HTAPS_MAX-2)
           || (unsigned)sy >= h - b_h - (HTAPS_MAX-2)){
            s->dsp.emulated_edge_mc(tmp + MB_SIZE, src, stride, b_w+HTAPS_MAX-1, b_h+HTAPS_MAX-1, sx, sy, w, h);
            src= tmp + MB_SIZE;
        }
//        assert(b_w == b_h || 2*b_w == b_h || b_w == 2*b_h);
//        assert(!(b_w&(b_w-1)));
        assert(b_w>1 && b_h>1);
        assert((tab_index>=0 && tab_index<4) || b_w==32);
        if((dx&3) || (dy&3) || !(b_w == b_h || 2*b_w == b_h || b_w == 2*b_h) || (b_w&(b_w-1)) || !s->plane[plane_index].fast_mc )
            mc_block(&s->plane[plane_index], dst, src, stride, b_w, b_h, dx, dy);
        else if(b_w==32){
            int y;
            for(y=0; y<b_h; y+=16){
                s->dsp.put_h264_qpel_pixels_tab[0][dy+(dx>>2)](dst + y*stride, src + 3 + (y+3)*stride,stride);
                s->dsp.put_h264_qpel_pixels_tab[0][dy+(dx>>2)](dst + 16 + y*stride, src + 19 + (y+3)*stride,stride);
            }
        }else if(b_w==b_h)
            s->dsp.put_h264_qpel_pixels_tab[tab_index  ][dy+(dx>>2)](dst,src + 3 + 3*stride,stride);
        else if(b_w==2*b_h){
            s->dsp.put_h264_qpel_pixels_tab[tab_index+1][dy+(dx>>2)](dst    ,src + 3       + 3*stride,stride);
            s->dsp.put_h264_qpel_pixels_tab[tab_index+1][dy+(dx>>2)](dst+b_h,src + 3 + b_h + 3*stride,stride);
        }else{
            assert(2*b_w==b_h);
            s->dsp.put_h264_qpel_pixels_tab[tab_index  ][dy+(dx>>2)](dst           ,src + 3 + 3*stride           ,stride);
            s->dsp.put_h264_qpel_pixels_tab[tab_index  ][dy+(dx>>2)](dst+b_w*stride,src + 3 + 3*stride+b_w*stride,stride);
        }
    }
}

#define mca(dx,dy,b_w)\
static void mc_block_hpel ## dx ## dy ## b_w(uint8_t *dst, const uint8_t *src, int stride, int h){\
    assert(h==b_w);\
    mc_block(NULL, dst, src-(HTAPS_MAX/2-1)-(HTAPS_MAX/2-1)*stride, stride, b_w, b_w, dx, dy);\
}

mca( 0, 0,16)
mca( 8, 0,16)
mca( 0, 8,16)
mca( 8, 8,16)
mca( 0, 0,8)
mca( 8, 0,8)
mca( 0, 8,8)
mca( 8, 8,8)

av_cold int ff_snow_common_init(AVCodecContext *avctx){
    SnowContext *s = avctx->priv_data;
    int width, height;
    int i, j;

    s->avctx= avctx;
    s->max_ref_frames=1; //just make sure its not an invalid value in case of no initial keyframe

    dsputil_init(&s->dsp, avctx);
    ff_dwt_init(&s->dwt);

#define mcf(dx,dy)\
    s->dsp.put_qpel_pixels_tab       [0][dy+dx/4]=\
    s->dsp.put_no_rnd_qpel_pixels_tab[0][dy+dx/4]=\
        s->dsp.put_h264_qpel_pixels_tab[0][dy+dx/4];\
    s->dsp.put_qpel_pixels_tab       [1][dy+dx/4]=\
    s->dsp.put_no_rnd_qpel_pixels_tab[1][dy+dx/4]=\
        s->dsp.put_h264_qpel_pixels_tab[1][dy+dx/4];

    mcf( 0, 0)
    mcf( 4, 0)
    mcf( 8, 0)
    mcf(12, 0)
    mcf( 0, 4)
    mcf( 4, 4)
    mcf( 8, 4)
    mcf(12, 4)
    mcf( 0, 8)
    mcf( 4, 8)
    mcf( 8, 8)
    mcf(12, 8)
    mcf( 0,12)
    mcf( 4,12)
    mcf( 8,12)
    mcf(12,12)

#define mcfh(dx,dy)\
    s->dsp.put_pixels_tab       [0][dy/4+dx/8]=\
    s->dsp.put_no_rnd_pixels_tab[0][dy/4+dx/8]=\
        mc_block_hpel ## dx ## dy ## 16;\
    s->dsp.put_pixels_tab       [1][dy/4+dx/8]=\
    s->dsp.put_no_rnd_pixels_tab[1][dy/4+dx/8]=\
        mc_block_hpel ## dx ## dy ## 8;

    mcfh(0, 0)
    mcfh(8, 0)
    mcfh(0, 8)
    mcfh(8, 8)

    init_qexp();

//    dec += FFMAX(s->chroma_h_shift, s->chroma_v_shift);

    width= s->avctx->width;
    height= s->avctx->height;

    s->spatial_idwt_buffer= av_mallocz(width*height*sizeof(IDWTELEM));
    s->spatial_dwt_buffer= av_mallocz(width*height*sizeof(DWTELEM)); //FIXME this does not belong here

    for(i=0; i<MAX_REF_FRAMES; i++)
        for(j=0; j<MAX_REF_FRAMES; j++)
            scale_mv_ref[i][j] = 256*(i+1)/(j+1);

    s->avctx->get_buffer(s->avctx, &s->mconly_picture);
    s->scratchbuf = av_malloc(s->mconly_picture.linesize[0]*7*MB_SIZE);

    return 0;
}

int ff_snow_common_init_after_header(AVCodecContext *avctx) {
    SnowContext *s = avctx->priv_data;
    int plane_index, level, orientation;

    for(plane_index=0; plane_index<3; plane_index++){
        int w= s->avctx->width;
        int h= s->avctx->height;

        if(plane_index){
            w>>= s->chroma_h_shift;
            h>>= s->chroma_v_shift;
        }
        s->plane[plane_index].width = w;
        s->plane[plane_index].height= h;

        for(level=s->spatial_decomposition_count-1; level>=0; level--){
            for(orientation=level ? 1 : 0; orientation<4; orientation++){
                SubBand *b= &s->plane[plane_index].band[level][orientation];

                b->buf= s->spatial_dwt_buffer;
                b->level= level;
                b->stride= s->plane[plane_index].width << (s->spatial_decomposition_count - level);
                b->width = (w + !(orientation&1))>>1;
                b->height= (h + !(orientation>1))>>1;

                b->stride_line = 1 << (s->spatial_decomposition_count - level);
                b->buf_x_offset = 0;
                b->buf_y_offset = 0;

                if(orientation&1){
                    b->buf += (w+1)>>1;
                    b->buf_x_offset = (w+1)>>1;
                }
                if(orientation>1){
                    b->buf += b->stride>>1;
                    b->buf_y_offset = b->stride_line >> 1;
                }
                b->ibuf= s->spatial_idwt_buffer + (b->buf - s->spatial_dwt_buffer);

                if(level)
                    b->parent= &s->plane[plane_index].band[level-1][orientation];
                //FIXME avoid this realloc
                av_freep(&b->x_coeff);
                b->x_coeff=av_mallocz(((b->width+1) * b->height+1)*sizeof(x_and_coeff));
            }
            w= (w+1)>>1;
            h= (h+1)>>1;
        }
    }

    return 0;
}

#define USE_HALFPEL_PLANE 0

static void halfpel_interpol(SnowContext *s, uint8_t *halfpel[4][4], AVFrame *frame){
    int p,x,y;

    for(p=0; p<3; p++){
        int is_chroma= !!p;
        int w= s->avctx->width  >>is_chroma;
        int h= s->avctx->height >>is_chroma;
        int ls= frame->linesize[p];
        uint8_t *src= frame->data[p];

        halfpel[1][p]= (uint8_t*)av_malloc(ls * (h+2*EDGE_WIDTH)) + EDGE_WIDTH*(1+ls);
        halfpel[2][p]= (uint8_t*)av_malloc(ls * (h+2*EDGE_WIDTH)) + EDGE_WIDTH*(1+ls);
        halfpel[3][p]= (uint8_t*)av_malloc(ls * (h+2*EDGE_WIDTH)) + EDGE_WIDTH*(1+ls);

        halfpel[0][p]= src;
        for(y=0; y<h; y++){
            for(x=0; x<w; x++){
                int i= y*ls + x;

                halfpel[1][p][i]= (20*(src[i] + src[i+1]) - 5*(src[i-1] + src[i+2]) + (src[i-2] + src[i+3]) + 16 )>>5;
            }
        }
        for(y=0; y<h; y++){
            for(x=0; x<w; x++){
                int i= y*ls + x;

                halfpel[2][p][i]= (20*(src[i] + src[i+ls]) - 5*(src[i-ls] + src[i+2*ls]) + (src[i-2*ls] + src[i+3*ls]) + 16 )>>5;
            }
        }
        src= halfpel[1][p];
        for(y=0; y<h; y++){
            for(x=0; x<w; x++){
                int i= y*ls + x;

                halfpel[3][p][i]= (20*(src[i] + src[i+ls]) - 5*(src[i-ls] + src[i+2*ls]) + (src[i-2*ls] + src[i+3*ls]) + 16 )>>5;
            }
        }

//FIXME border!
    }
}

void ff_snow_release_buffer(AVCodecContext *avctx)
{
    SnowContext *s = avctx->priv_data;
    int i;

    if(s->last_picture[s->max_ref_frames-1].data[0]){
        avctx->release_buffer(avctx, &s->last_picture[s->max_ref_frames-1]);
        for(i=0; i<9; i++)
            if(s->halfpel_plane[s->max_ref_frames-1][1+i/3][i%3])
                av_free(s->halfpel_plane[s->max_ref_frames-1][1+i/3][i%3] - EDGE_WIDTH*(1+s->current_picture.linesize[i%3]));
    }
}

int ff_snow_frame_start(SnowContext *s){
   AVFrame tmp;
   int w= s->avctx->width; //FIXME round up to x16 ?
   int h= s->avctx->height;

    if (s->current_picture.data[0] && !(s->avctx->flags&CODEC_FLAG_EMU_EDGE)) {
        s->dsp.draw_edges(s->current_picture.data[0],
                          s->current_picture.linesize[0], w   , h   ,
                          EDGE_WIDTH  , EDGE_WIDTH  , EDGE_TOP | EDGE_BOTTOM);
        s->dsp.draw_edges(s->current_picture.data[1],
                          s->current_picture.linesize[1], w>>1, h>>1,
                          EDGE_WIDTH/2, EDGE_WIDTH/2, EDGE_TOP | EDGE_BOTTOM);
        s->dsp.draw_edges(s->current_picture.data[2],
                          s->current_picture.linesize[2], w>>1, h>>1,
                          EDGE_WIDTH/2, EDGE_WIDTH/2, EDGE_TOP | EDGE_BOTTOM);
    }

    ff_snow_release_buffer(s->avctx);

    tmp= s->last_picture[s->max_ref_frames-1];
    memmove(s->last_picture+1, s->last_picture, (s->max_ref_frames-1)*sizeof(AVFrame));
    memmove(s->halfpel_plane+1, s->halfpel_plane, (s->max_ref_frames-1)*sizeof(void*)*4*4);
    if(USE_HALFPEL_PLANE && s->current_picture.data[0])
        halfpel_interpol(s, s->halfpel_plane[0], &s->current_picture);
    s->last_picture[0]= s->current_picture;
    s->current_picture= tmp;

    if(s->keyframe){
        s->ref_frames= 0;
    }else{
        int i;
        for(i=0; i<s->max_ref_frames && s->last_picture[i].data[0]; i++)
            if(i && s->last_picture[i-1].key_frame)
                break;
        s->ref_frames= i;
        if(s->ref_frames==0){
            av_log(s->avctx,AV_LOG_ERROR, "No reference frames\n");
            return -1;
        }
    }

    s->current_picture.reference= 1;
    if(s->avctx->get_buffer(s->avctx, &s->current_picture) < 0){
        av_log(s->avctx, AV_LOG_ERROR, "get_buffer() failed\n");
        return -1;
    }

    s->current_picture.key_frame= s->keyframe;

    return 0;
}

av_cold void ff_snow_common_end(SnowContext *s)
{
    int plane_index, level, orientation, i;

    av_freep(&s->spatial_dwt_buffer);
    av_freep(&s->spatial_idwt_buffer);

    s->m.me.temp= NULL;
    av_freep(&s->m.me.scratchpad);
    av_freep(&s->m.me.map);
    av_freep(&s->m.me.score_map);
    av_freep(&s->m.obmc_scratchpad);

    av_freep(&s->block);
    av_freep(&s->scratchbuf);

    for(i=0; i<MAX_REF_FRAMES; i++){
        av_freep(&s->ref_mvs[i]);
        av_freep(&s->ref_scores[i]);
        if(s->last_picture[i].data[0])
            s->avctx->release_buffer(s->avctx, &s->last_picture[i]);
    }

    for(plane_index=0; plane_index<3; plane_index++){
        for(level=s->spatial_decomposition_count-1; level>=0; level--){
            for(orientation=level ? 1 : 0; orientation<4; orientation++){
                SubBand *b= &s->plane[plane_index].band[level][orientation];

                av_freep(&b->x_coeff);
            }
        }
    }
    if (s->mconly_picture.data[0])
        s->avctx->release_buffer(s->avctx, &s->mconly_picture);
    if (s->current_picture.data[0])
        s->avctx->release_buffer(s->avctx, &s->current_picture);
}

