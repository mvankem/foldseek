/*
=============================================================
   Implementation of TM-align in C/C++

   This program is written by Jianyi Yang at
   Yang Zhang lab
   And it is updated by Jianjie Wu at
   Yang Zhang lab
   Department of Computational Medicine and Bioinformatics
   University of Michigan
   100 Washtenaw Avenue, Ann Arbor, MI 48109-2218

   Please report bugs and questions to zhng@umich.edu
=============================================================
*/

#include "NW.h"
#include "Kabsch.h"
#include "affineneedlemanwunsch.h"
#include "basic_fun.h"



void parameter_set4search(const int xlen, const int ylen,
                          float &D0_MIN, float &Lnorm,
                          float &score_d8, float &d0, float &d0_search, float &dcu0)
{
    //parameter initilization for searching: D0_MIN, Lnorm, d0, d0_search, score_d8
    D0_MIN=0.5;
    dcu0=4.25;                       //update 3.85-->4.25

    Lnorm=std::min(xlen, ylen);        //normaliz TMscore by this in searching
    if (Lnorm<=19){                    //update 15-->19
        d0=0.168;                   //update 0.5-->0.168
    } else {
        d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    }
    D0_MIN=d0+0.8;              //this should be moved to above
    d0=D0_MIN;                  //update: best for search

    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;

    score_d8=1.5*pow(Lnorm*1.0, 0.3)+3.5; //remove pairs with dis>d8 during search & final
}

void parameter_set4final(const float len, float &D0_MIN, float &Lnorm,
                         float &d0, float &d0_search)
{
    D0_MIN=0.5;

    Lnorm=len;            //normaliz TMscore by this in searching
    if (Lnorm<=21) d0=0.5;
    else d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    if (d0<D0_MIN) d0=D0_MIN;
    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;
}

void parameter_set4scale(const int len, const float d_s, float &Lnorm,
                         float &d0, float &d0_search)
{
    d0=d_s;
    Lnorm=len;            //normaliz TMscore by this in searching
    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;
}

//     1, collect those residues with dis<d;
//     2, calculate TMscore
//int score_fun8( Coordinates &xa, Coordinates &ya, int n_ali, float d, int i_ali[],
//                float *score1, int score_sum_method, const float Lnorm,
//                const float score_d8, const float d0
//                , float * mem
//
//                )
//{
//    float score_sum=0, di;
//    float d_tmp=d*d;
//    float d02=d0*d0;
//    float score_d8_cut = score_d8*score_d8;
//
//    int i, n_cut, inc=0;
//    float *BasicFunction::distArray = mem;
//    while(1)
//    {
//        n_cut=0;
//        score_sum=0;
//        simd_float sum = simdf32_set(0);
//        simd_float vscore_d8_cut = simdf32_set(score_d8_cut);
//        simd_float vd02 = simdf32_set(d02);
//        simd_float one = simdf32_set(1.0f);
//        for(i=0; i < n_ali; i+=VECSIZE_FLOAT){
//            //    float d1=xx-yx;
//            //    float d2=xy-yy;
//            //    float d3=xz-yz;
//            //    return (d1*d1 + d2*d2 + d3*d3);
//            simd_float xa_x = simdf32_load(&xa.x[i]);
//            simd_float ya_x = simdf32_load(&ya.x[i]);
//            simd_float xa_y = simdf32_load(&xa.y[i]);
//            simd_float ya_y = simdf32_load(&ya.y[i]);
//            simd_float xa_z = simdf32_load(&xa.z[i]);
//            simd_float ya_z = simdf32_load(&ya.z[i]);
//            ya_x = simdf32_sub(xa_x, ya_x);
//            ya_y = simdf32_sub(xa_y, ya_y);
//            ya_z = simdf32_sub(xa_z, ya_z);
//            ya_x = simdf32_mul(ya_x, ya_x);
//            ya_y = simdf32_mul(ya_y, ya_y);
//            ya_z = simdf32_mul(ya_z, ya_z);
//            simd_float res = simdf32_add(ya_x, ya_y);
//            simd_float di = simdf32_add(res, ya_z);
//            simdf32_store(&BasicFunction::distArray[i], di);
//            simd_float di_lt_score_d8 = simdf32_lt(di, vscore_d8_cut);
//            simd_float oneDividedDist = simdf32_div(one, simdf32_add(one, simdf32_div(di,vd02)));
//            sum = simdf32_add(sum, (simd_float)simdi_and((simd_int) di_lt_score_d8, (simd_int) oneDividedDist ));
//        }
//        for(i=0; i < VECSIZE_FLOAT; i++){
//            score_sum+=((float*)&sum)[i];
//        }
//
//        for(i=0; i<n_ali; i++)
//        {
//            di = BasicFunction::distArray[i];
//            i_ali[n_cut]=i;
//            n_cut+=(di<d_tmp);
//            //score_sum += (di<=score_d8_cut) ? 1/(1+di/d02) : 0;
//            //else score_sum += 1/(1+di/d02);
//        }
//        //there are not enough feasible pairs, reliefe the threshold
//        if(n_cut<3 && n_ali>3)
//        {
//            inc++;
//            double dinc=(d+inc*0.5);
//            d_tmp = dinc * dinc;
//        }
//        else break;
//    }
//
//    *score1=score_sum/Lnorm;
//    return n_cut;
//}


//     1, collect those residues with dis<d;
//     2, calculate TMscore
int score_fun8( Coordinates &xa, Coordinates &ya, int n_ali, float d, int i_ali[],
                float *score1, const float Lnorm,
                const float score_d8, const float d0, float * mem)
{
    float score_sum=0, di;
    float d_tmp=d*d;
    float d02=d0*d0;
    float score_d8_cut = score_d8*score_d8;
    int i, n_cut, inc=0;
    float *distArray = mem;
    float *sumArray = mem+((n_ali/VECSIZE_FLOAT+1)*VECSIZE_FLOAT);

    while(1)
    {
        n_cut=0;
        score_sum=0;
        simd_float vscore_d8_cut = simdf32_set(score_d8_cut);
        simd_float vd02 = simdf32_set(d02);
        simd_float one = simdf32_set(1.0f);
        for(i=0; i < n_ali; i+=VECSIZE_FLOAT){
            //    float d1=xx-yx;
            //    float d2=xy-yy;
            //    float d3=xz-yz;
            //    return (d1*d1 + d2*d2 + d3*d3);
            simd_float xa_x = simdf32_load(&xa.x[i]);
            simd_float ya_x = simdf32_load(&ya.x[i]);
            simd_float xa_y = simdf32_load(&xa.y[i]);
            simd_float ya_y = simdf32_load(&ya.y[i]);
            simd_float xa_z = simdf32_load(&xa.z[i]);
            simd_float ya_z = simdf32_load(&ya.z[i]);
            ya_x = simdf32_sub(xa_x, ya_x);
            ya_y = simdf32_sub(xa_y, ya_y);
            ya_z = simdf32_sub(xa_z, ya_z);
            ya_x = simdf32_mul(ya_x, ya_x);
            ya_y = simdf32_mul(ya_y, ya_y);
            ya_z = simdf32_mul(ya_z, ya_z);
            simd_float res = simdf32_add(ya_x, ya_y);
            simd_float di = simdf32_add(res, ya_z);
            simdf32_store(&distArray[i], di);
            simd_float di_lt_score_d8 = simdf32_lt(di, vscore_d8_cut);
            simd_float oneDividedDist = simdf32_div(one, simdf32_add(one, simdf32_div(di,vd02)));
            //sum = simdf32_add(sum, (simd_float));
            simdf32_store(&sumArray[i], (simd_float)simdi_and((simd_int) di_lt_score_d8, (simd_int) oneDividedDist ));
        }


        for(i=0; i<n_ali; i++)
        {
            di = distArray[i];
            i_ali[n_cut]=i;
            n_cut+=(di<d_tmp);
            score_sum+=sumArray[i];
            //score_sum += (di<=score_d8_cut) ? 1/(1+di/d02) : 0;
            //else score_sum += 1/(1+di/d02);
        }
        //there are not enough feasible pairs, reliefe the threshold
        if(n_cut<3 && n_ali>3)
        {
            inc++;
            double dinc=(d+inc*0.5);
            d_tmp = dinc * dinc;
        }
        else break;
    }

    *score1=score_sum/Lnorm;
    return n_cut;
}


//int score_fun8( Coordinates &xa, Coordinates &ya, int n_ali, float d, int i_ali[],
//                float *score1, int score_sum_method, const float Lnorm,
//                const float score_d8, const float d0, float * mem)
//{
//    float score_sum=0, di;
//    float d_tmp=d*d;
//    float d02=d0*d0;
//    float score_d8_cut = score_d8*score_d8;
//
//    int i, n_cut, inc=0;
//
//    while(1)
//    {
//        n_cut=0;
//        score_sum=0;
//        for(i=0; i<n_ali; i++)
//        {
//            di = BasicFunction::dist(xa.x[i], xa.y[i], xa.z[i], ya.x[i], ya.y[i], ya.z[i]);
//            if(di<d_tmp)
//            {
//                i_ali[n_cut]=i;
//                n_cut++;
//            }
//            if(score_sum_method==8)
//            {
//                if(di<=score_d8_cut) score_sum += 1/(1+di/d02);
//            }
//            else score_sum += 1/(1+di/d02);
//        }
//        //there are not enough feasible pairs, reliefe the threshold
//        if(n_cut<3 && n_ali>3)
//        {
//            inc++;
//            double dinc=(d+inc*0.5);
//            d_tmp = dinc * dinc;
//        }
//        else break;
//    }
//
//    *score1=score_sum/Lnorm;
//    return n_cut;
//}

int score_fun8_standard(Coordinates &xa, Coordinates &ya, int n_ali, float d,
                        int i_ali[], float *score1, int score_sum_method,
                        float score_d8, float d0)
{
    float score_sum = 0, di;
    float d_tmp = d*d;
    float d02 = d0*d0;
    float score_d8_cut = score_d8*score_d8;

    int i, n_cut, inc = 0;
    while (1)
    {
        n_cut = 0;
        score_sum = 0;
        for (i = 0; i<n_ali; i++)
        {
            di = BasicFunction::BasicFunction::dist(xa.x[i], xa.y[i], xa.z[i],
                      ya.x[i], ya.y[i], ya.z[i]);
            if (di<d_tmp)
            {
                i_ali[n_cut] = i;
                n_cut++;
            }
            if (score_sum_method == 8)
            {
                if (di <= score_d8_cut) score_sum += 1 / (1 + di / d02);
            }
            else
            {
                score_sum += 1 / (1 + di / d02);
            }
        }
        //there are not enough feasible pairs, reliefe the threshold
        if (n_cut<3 && n_ali>3)
        {
            inc++;
            double dinc = (d + inc*0.5);
            d_tmp = dinc * dinc;
        }
        else break;
    }

    *score1 = score_sum / n_ali;
    return n_cut;
}

/*bool Kabsch(Coordinates & x, Coordinates & y, int n, int mode, double *rms,
            double t[3], double u[3][3])
{
    int i, j, m, m1, l, k;
    double e0, rms1, d, h, g;
    double cth, sth, sqrth, p, det, sigma;
    double xc[3], yc[3];
    double a[3][3], b[3][3], r[3][3], e[3], rr[6], ss[6];
    double sqrt3 = 1.73205080756888, tol = 0.01;
    int ip[] = { 0, 1, 3, 1, 2, 4, 3, 4, 5 };
    int ip2312[] = { 1, 2, 0, 1 };

    int a_failed = 0, b_failed = 0;
    double epsilon = 0.00000001;

    //initialization
    *rms = 0;
    rms1 = 0;
    e0 = 0;
    double c1[3], c2[3];
    double s1[3], s2[3];
    double sx[3], sy[3], sz[3];
    for (i = 0; i < 3; i++)
    {
        s1[i] = 0.0;
        s2[i] = 0.0;

        sx[i] = 0.0;
        sy[i] = 0.0;
        sz[i] = 0.0;
    }

    for (i = 0; i<3; i++)
    {
        xc[i] = 0.0;
        yc[i] = 0.0;
        t[i] = 0.0;
        for (j = 0; j<3; j++)
        {
            u[i][j] = 0.0;
            r[i][j] = 0.0;
            a[i][j] = 0.0;
            if (i == j)
            {
                u[i][j] = 1.0;
                a[i][j] = 1.0;
            }
        }
    }

    if (n<1) return false;

    //compute centers for vector sets x, y
    for (i = 0; i<n; i++)
    {
        for (j = 0; j < 3; j++)
        {
            c1[0] = x.x[i];
            c1[1] = x.y[i];
            c1[2] = x.z[i];

            c2[0] = y.x[i];
            c2[1] = y.y[i];
            c2[2] = y.z[i];

            s1[j] += c1[0]+c1[1]+c1[2];
            s2[j] += c2[0]+c2[1]+c2[2];
        }

        for (j = 0; j < 3; j++)
        {
            sx[j] += c1[0] * c2[j];
            sy[j] += c1[1] * c2[j];
            sz[j] += c1[2] * c2[j];
        }
    }
    for (i = 0; i < 3; i++)
    {
        xc[i] = s1[i] / n;
        yc[i] = s2[i] / n;
    }
    if (mode == 2 || mode == 0)
        for (int mm = 0; mm < n; mm++){
            e0 += (x.x[mm] - xc[0]) * (x.x[mm] - xc[0]) +
                  (y.x[mm] - yc[0]) * (y.x[mm] - yc[0]);
            e0 += (x.y[mm] - xc[1]) * (x.y[mm] - xc[1]) +
                  (y.y[mm] - yc[1]) * (y.y[mm] - yc[1]);
            e0 += (x.z[mm] - xc[2]) * (x.z[mm] - xc[2]) +
                  (y.z[mm] - yc[2]) * (y.z[mm] - yc[2]);
        }
    for (j = 0; j < 3; j++)
    {
        r[j][0] = sx[j] - s1[0] * s2[j] / n;
        r[j][1] = sy[j] - s1[1] * s2[j] / n;
        r[j][2] = sz[j] - s1[2] * s2[j] / n;
    }

    //compute determinant of matrix r
    det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])\
        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])\
        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
    sigma = det;

    //compute tras(r)*r
    m = 0;
    for (j = 0; j<3; j++)
    {
        for (i = 0; i <= j; i++)
        {
            rr[m] = r[0][i] * r[0][j] + r[1][i] * r[1][j] + r[2][i] * r[2][j];
            m++;
        }
    }

    double spur = (rr[0] + rr[2] + rr[5]) / 3.0;
    double cof = (((((rr[2] * rr[5] - rr[4] * rr[4]) + rr[0] * rr[5])\
        - rr[3] * rr[3]) + rr[0] * rr[2]) - rr[1] * rr[1]) / 3.0;
    det = det*det;

    for (i = 0; i<3; i++) e[i] = spur;

    if (spur>0)
    {
        d = spur*spur;
        h = d - cof;
        g = (spur*cof - det) / 2.0 - spur*h;

        if (h>0)
        {
            sqrth = sqrt(h);
            d = h*h*h - g*g;
            if (d<0.0) d = 0.0;
            d = atan2(sqrt(d), -g) / 3.0;
            cth = sqrth * cos(d);
            sth = sqrth*sqrt3*sin(d);
            e[0] = (spur + cth) + cth;
            e[1] = (spur - cth) + sth;
            e[2] = (spur - cth) - sth;

            if (mode != 0)
            {//compute a
                for (l = 0; l<3; l = l + 2)
                {
                    d = e[l];
                    ss[0] = (d - rr[2]) * (d - rr[5]) - rr[4] * rr[4];
                    ss[1] = (d - rr[5]) * rr[1] + rr[3] * rr[4];
                    ss[2] = (d - rr[0]) * (d - rr[5]) - rr[3] * rr[3];
                    ss[3] = (d - rr[2]) * rr[3] + rr[1] * rr[4];
                    ss[4] = (d - rr[0]) * rr[4] + rr[1] * rr[3];
                    ss[5] = (d - rr[0]) * (d - rr[2]) - rr[1] * rr[1];

                    if (fabs(ss[0]) <= epsilon) ss[0] = 0.0;
                    if (fabs(ss[1]) <= epsilon) ss[1] = 0.0;
                    if (fabs(ss[2]) <= epsilon) ss[2] = 0.0;
                    if (fabs(ss[3]) <= epsilon) ss[3] = 0.0;
                    if (fabs(ss[4]) <= epsilon) ss[4] = 0.0;
                    if (fabs(ss[5]) <= epsilon) ss[5] = 0.0;

                    if (fabs(ss[0]) >= fabs(ss[2]))
                    {
                        j = 0;
                        if (fabs(ss[0]) < fabs(ss[5])) j = 2;
                    }
                    else if (fabs(ss[2]) >= fabs(ss[5])) j = 1;
                    else j = 2;

                    d = 0.0;
                    j = 3 * j;
                    for (i = 0; i<3; i++)
                    {
                        k = ip[i + j];
                        a[i][l] = ss[k];
                        d = d + ss[k] * ss[k];
                    }


                    //if( d > 0.0 ) d = 1.0 / sqrt(d);
                    if (d > epsilon) d = 1.0 / sqrt(d);
                    else d = 0.0;
                    for (i = 0; i<3; i++) a[i][l] = a[i][l] * d;
                }//for l

                d = a[0][0] * a[0][2] + a[1][0] * a[1][2] + a[2][0] * a[2][2];
                if ((e[0] - e[1]) >(e[1] - e[2]))
                {
                    m1 = 2;
                    m = 0;
                }
                else
                {
                    m1 = 0;
                    m = 2;
                }
                p = 0;
                for (i = 0; i<3; i++)
                {
                    a[i][m1] = a[i][m1] - d*a[i][m];
                    p = p + a[i][m1] * a[i][m1];
                }
                if (p <= tol)
                {
                    p = 1.0;
                    for (i = 0; i<3; i++)
                    {
                        if (p < fabs(a[i][m])) continue;
                        p = fabs(a[i][m]);
                        j = i;
                    }
                    k = ip2312[j];
                    l = ip2312[j + 1];
                    p = sqrt(a[k][m] * a[k][m] + a[l][m] * a[l][m]);
                    if (p > tol)
                    {
                        a[j][m1] = 0.0;
                        a[k][m1] = -a[l][m] / p;
                        a[l][m1] = a[k][m] / p;
                    }
                    else a_failed = 1;
                }//if p<=tol
                else
                {
                    p = 1.0 / sqrt(p);
                    for (i = 0; i<3; i++) a[i][m1] = a[i][m1] * p;
                }//else p<=tol
                if (a_failed != 1)
                {
                    a[0][1] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
                    a[1][1] = a[2][2] * a[0][0] - a[2][0] * a[0][2];
                    a[2][1] = a[0][2] * a[1][0] - a[0][0] * a[1][2];
                }
            }//if(mode!=0)
        }//h>0

        //compute b anyway
        if (mode != 0 && a_failed != 1)//a is computed correctly
        {
            //compute b
            for (l = 0; l<2; l++)
            {
                d = 0.0;
                for (i = 0; i<3; i++)
                {
                    b[i][l] = r[i][0] * a[0][l] +
                              r[i][1] * a[1][l] + r[i][2] * a[2][l];
                    d = d + b[i][l] * b[i][l];
                }
                //if( d > 0 ) d = 1.0 / sqrt(d);
                if (d > epsilon) d = 1.0 / sqrt(d);
                else d = 0.0;
                for (i = 0; i<3; i++) b[i][l] = b[i][l] * d;
            }
            d = b[0][0] * b[0][1] + b[1][0] * b[1][1] + b[2][0] * b[2][1];
            p = 0.0;

            for (i = 0; i<3; i++)
            {
                b[i][1] = b[i][1] - d*b[i][0];
                p += b[i][1] * b[i][1];
            }

            if (p <= tol)
            {
                p = 1.0;
                for (i = 0; i<3; i++)
                {
                    if (p<fabs(b[i][0])) continue;
                    p = fabs(b[i][0]);
                    j = i;
                }
                k = ip2312[j];
                l = ip2312[j + 1];
                p = sqrt(b[k][0] * b[k][0] + b[l][0] * b[l][0]);
                if (p > tol)
                {
                    b[j][1] = 0.0;
                    b[k][1] = -b[l][0] / p;
                    b[l][1] = b[k][0] / p;
                }
                else b_failed = 1;
            }//if( p <= tol )
            else
            {
                p = 1.0 / sqrt(p);
                for (i = 0; i<3; i++) b[i][1] = b[i][1] * p;
            }
            if (b_failed != 1)
            {
                b[0][2] = b[1][0] * b[2][1] - b[1][1] * b[2][0];
                b[1][2] = b[2][0] * b[0][1] - b[2][1] * b[0][0];
                b[2][2] = b[0][0] * b[1][1] - b[0][1] * b[1][0];
                //compute u
                for (i = 0; i<3; i++)
                    for (j = 0; j<3; j++)
                        u[i][j] = b[i][0] * a[j][0] +
                                  b[i][1] * a[j][1] + b[i][2] * a[j][2];
            }

            //compute t
            for (i = 0; i<3; i++)
                t[i] = ((yc[i] - u[i][0] * xc[0]) - u[i][1] * xc[1]) -
                       u[i][2] * xc[2];
        }//if(mode!=0 && a_failed!=1)
    }//spur>0
    else //just compute t and errors
    {
        //compute t
        for (i = 0; i<3; i++)
            t[i] = ((yc[i] - u[i][0] * xc[0]) - u[i][1] * xc[1]) -
                   u[i][2] * xc[2];
    }//else spur>0

    //compute rms
    for (i = 0; i<3; i++)
    {
        if (e[i] < 0) e[i] = 0;
        e[i] = sqrt(e[i]);
    }
    d = e[2];
    if (sigma < 0.0) d = -d;
    d = (d + e[1]) + e[0];

    if (mode == 2 || mode == 0)
    {
        rms1 = (e0 - d) - d;
        if (rms1 < 0.0) rms1 = 0.0;
    }

    *rms = rms1;
    return true;
}
*/
bool KabschFast(Coordinates & x,
                 Coordinates & y,
                 int n,
                 float *rms,
                 float t[3],
                 float u[3][3],
                 float * mem){
    float r[16] __attribute__ ((aligned (16)));

//    Kabsch(x, y, n, 2, rms,
//            t, u);
    *rms = kabsch_quat_soa_avx(n, x.x, x.y, x.z, y.x, y.y, y.z, r, mem);
    t[0] = r[3];
    t[1] = r[7];
    t[2] = r[11];
    u[0][0] = r[0];
    u[0][1] = r[1];
    u[0][2] = r[2];
    u[1][0] = r[4];
    u[1][1] = r[5];
    u[1][2] = r[6];
    u[2][0] = r[8];
    u[2][1] = r[9];
    u[2][2] = r[10];
    return true;
}

double TMscore8_search(Coordinates &r1, Coordinates &r2,
                       Coordinates &xtm, Coordinates & ytm,
                       Coordinates &xt, int Lali, float t0[3], float u0[3][3], int simplify_step,
                       float *Rcomm, float local_d0_search, float Lnorm,
                       float score_d8, float d0, float * mem)
{
    int i, m;
    float score_max, score, rmsd;
    const int kmax=Lali;
    int k_ali[kmax], ka, k;
    float t[3];
    float u[3][3];
    float d;


    //iterative parameters
    int n_it=10;            //maximum number of iterations
    int n_init_max=6; //maximum number of different fragment length
    int L_ini[n_init_max];  //fragment lengths, Lali, Lali/2, Lali/4 ... 4
    int L_ini_min=4;
    if(Lali<L_ini_min) L_ini_min=Lali;

    int n_init=0, i_init;
    for(i=0; i<n_init_max-1; i++)
    {
        n_init++;
        L_ini[i]=(int) (Lali/pow(2.0, (double) i));
        if(L_ini[i]<=L_ini_min)
        {
            L_ini[i]=L_ini_min;
            break;
        }
    }
    if(i==n_init_max-1)
    {
        n_init++;
        L_ini[i]=L_ini_min;
    }

    score_max=-1;
    //find the maximum score starting from local structures superposition
    int i_ali[kmax], n_cut;
    int L_frag; //fragment length
    int iL_max; //maximum starting postion for the fragment

    for(i_init=0; i_init<n_init; i_init++)
    {
        L_frag=L_ini[i_init];
        iL_max=Lali-L_frag;

        i=0;
        while(1)
        {
            //extract the fragment starting from position i
            ka=0;
            for(k=0; k<L_frag; k++)
            {
                int kk=k+i;
                r1.x[k]=xtm.x[kk];
                r1.y[k]=xtm.y[kk];
                r1.z[k]=xtm.z[kk];

                r2.x[k]=ytm.x[kk];
                r2.y[k]=ytm.y[kk];
                r2.z[k]=ytm.z[kk];

                k_ali[ka]=kk;
                ka++;
            }

            //extract rotation matrix based on the fragment
            //float r[16] __attribute__ ((aligned (16)));
            //float rmsdbla = kabsch_quat_soa_sse2(L_frag, NULL, r1.x, r1.y, r1.z, r2.x, r2.y, r2.z, r);
            //float rmsdbla;
//            float bR[3][3];
//            float bt[3];
            //float TMscore = tmscore_cpu_soa_sse2(L_frag, r1.x, r1.y, r1.z, r2.x, r2.y, r2.z, bR, bt, &rmsdbla);
            KabschFast(r1, r2, L_frag, &rmsd, t, u, mem);
            //Kabsch(r1, r2, L_frag, 1, &rmsd, t, u);//, mem);

            if (simplify_step != 1)
                *Rcomm = 0;
            BasicFunction::do_rotation(xtm, xt, Lali, t, u);

            //get subsegment of this fragment
            d = local_d0_search - 1;
            n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score,
                             Lnorm, score_d8, d0, mem);
            if(score>score_max)
            {
                score_max=score;

                //save the rotation matrix
                for(k=0; k<3; k++)
                {
                    t0[k]=t[k];
                    u0[k][0]=u[k][0];
                    u0[k][1]=u[k][1];
                    u0[k][2]=u[k][2];
                }
            }

            //try to extend the alignment iteratively
            d = local_d0_search + 1;
            for(int it=0; it<n_it; it++)
            {
                ka=0;
                for(k=0; k<n_cut; k++)
                {
                    m=i_ali[k];
                    r1.x[k]=xtm.x[m];
                    r1.y[k]=xtm.y[m];
                    r1.z[k]=xtm.z[m];

                    r2.x[k]=ytm.x[m];
                    r2.y[k]=ytm.y[m];
                    r2.z[k]=ytm.z[m];

                    k_ali[ka]=m;
                    ka++;
                }
                //KabschFast();
                //extract rotation matrix based on the fragment
                //float rmsdbla = kabsch_quat_soa_sse2(L_frag, NULL, r1.x, r1.y, r1.z, r2.x, r2.y, r2.z, r);

                KabschFast(r1, r2, n_cut, &rmsd, t, u, mem);
//                Kabsch(r1, r2, n_cut, 1, &rmsd, t, u);

                BasicFunction::do_rotation(xtm, xt, Lali, t, u);
                n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score,
                                 Lnorm, score_d8, d0, mem);
                if(score>score_max)
                {
                    score_max=score;

                    //save the rotation matrix
                    for(k=0; k<3; k++)
                    {
                        t0[k]=t[k];
                        u0[k][0]=u[k][0];
                        u0[k][1]=u[k][1];
                        u0[k][2]=u[k][2];
                    }
                }

                //check if it converges
                if(n_cut==ka)
                {
                    for(k=0; k<n_cut; k++)
                    {
                        if(i_ali[k]!=k_ali[k]) break;
                    }
                    if(k==n_cut) break;
                }
            } //for iteration

            if(i<iL_max)
            {
                i=i+simplify_step; //shift the fragment
                if(i>iL_max) i=iL_max;  //do this to use the last missed fragment
            }
            else if(i>=iL_max) break;
        }//while(1)
        //end of one fragment
    }//for(i_init
    return score_max;
}


double TMscore8_search_standard(Coordinates &r1, Coordinates &r2,
                                Coordinates &xtm, Coordinates &ytm, Coordinates &xt, int Lali,
                                float t0[3], float u0[3][3], int simplify_step, int score_sum_method,
                                float *Rcomm, float local_d0_search, float score_d8, float d0, float * mem)
{
    int i, m;
    float score_max, score, rmsd;
    const int kmax = Lali;
    int k_ali[kmax], ka, k;
    float t[3];
    float u[3][3];
    float d;

    //iterative parameters
    int n_it = 20;            //maximum number of iterations
    int n_init_max = 6; //maximum number of different fragment length
    int L_ini[n_init_max];  //fragment lengths, Lali, Lali/2, Lali/4 ... 4
    int L_ini_min = 4;
    if (Lali<L_ini_min) L_ini_min = Lali;

    int n_init = 0, i_init;
    for (i = 0; i<n_init_max - 1; i++)
    {
        n_init++;
        L_ini[i] = (int)(Lali / pow(2.0, (double)i));
        if (L_ini[i] <= L_ini_min)
        {
            L_ini[i] = L_ini_min;
            break;
        }
    }
    if (i == n_init_max - 1)
    {
        n_init++;
        L_ini[i] = L_ini_min;
    }

    score_max = -1;
    //find the maximum score starting from local structures superposition
    int i_ali[kmax], n_cut;
    int L_frag; //fragment length
    int iL_max; //maximum starting postion for the fragment

    for (i_init = 0; i_init<n_init; i_init++)
    {
        L_frag = L_ini[i_init];
        iL_max = Lali - L_frag;

        i = 0;
        while (1)
        {
            //extract the fragment starting from position i
            ka = 0;
            for (k = 0; k<L_frag; k++)
            {
                int kk = k + i;
                r1.x[k] = xtm.x[kk];
                r1.y[k] = xtm.y[kk];
                r1.z[k] = xtm.z[kk];

                r2.x[k] = ytm.x[kk];
                r2.y[k] = ytm.y[kk];
                r2.z[k] = ytm.z[kk];

                k_ali[ka] = kk;
                ka++;
            }
            //extract rotation matrix based on the fragment
            KabschFast(r1, r2, L_frag, &rmsd, t, u, mem);
            if (simplify_step != 1)
                *Rcomm = 0;
            BasicFunction::do_rotation(xtm, xt, Lali, t, u);

            //get subsegment of this fragment
            d = local_d0_search - 1;
            n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, &score,
                                        score_sum_method, score_d8, d0);

            if (score>score_max)
            {
                score_max = score;

                //save the rotation matrix
                for (k = 0; k<3; k++)
                {
                    t0[k] = t[k];
                    u0[k][0] = u[k][0];
                    u0[k][1] = u[k][1];
                    u0[k][2] = u[k][2];
                }
            }

            //try to extend the alignment iteratively
            d = local_d0_search + 1;
            for (int it = 0; it<n_it; it++)
            {
                ka = 0;
                for (k = 0; k<n_cut; k++)
                {
                    m = i_ali[k];
                    r1.x[k] = xtm.x[m];
                    r1.y[k] = xtm.y[m];
                    r1.z[k] = xtm.z[m];

                    r2.x[k] = ytm.x[m];
                    r2.y[k] = ytm.y[m];
                    r2.z[k] = ytm.z[m];

                    k_ali[ka] = m;
                    ka++;
                }
                //extract rotation matrix based on the fragment
                KabschFast(r1, r2, n_cut, &rmsd, t, u, mem);
                BasicFunction::do_rotation(xtm, xt, Lali, t, u);
                n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, &score,
                                            score_sum_method, score_d8, d0);
                if (score>score_max)
                {
                    score_max = score;

                    //save the rotation matrix
                    for (k = 0; k<3; k++)
                    {
                        t0[k] = t[k];
                        u0[k][0] = u[k][0];
                        u0[k][1] = u[k][1];
                        u0[k][2] = u[k][2];
                    }
                }

                //check if it converges
                if (n_cut == ka)
                {
                    for (k = 0; k<n_cut; k++)
                    {
                        if (i_ali[k] != k_ali[k]) break;
                    }
                    if (k == n_cut) break;
                }
            } //for iteration

            if (i<iL_max)
            {
                i = i + simplify_step; //shift the fragment
                if (i>iL_max) i = iL_max;  //do this to use the last missed fragment
            }
            else if (i >= iL_max) break;
        }//while(1)
        //end of one fragment
    }//for(i_init
    return score_max;
}


double get_score4pareun(Coordinates &r1, Coordinates &r2, Coordinates &xtm, Coordinates &ytm,
                        const Coordinates &x, const Coordinates &y, int * queryToTargetMapping,
                        int queryLen,
                        float t[3], float u[3][3], float * mem ) {
    float rms, tmscore;
    int i, j, k;
    float d0;

    float D0_MIN = 0.5;
    int Lnorm = queryLen;

    if (Lnorm > 21)
        d0 = (1.24*pow((Lnorm*1.0 - 15), 1.0 / 3) - 1.8);
    else
        d0 = D0_MIN;
    if (d0 < D0_MIN)
        d0 = D0_MIN;

    k=0;
    for(j=0; j< queryLen; j++)
    {
        i=queryToTargetMapping[j];
        //i = invmap[k][0];
        //j = invmap[k][1];
        if(i>=0)
        {
            r1.x[k]=x.x[i];
            r1.y[k]=x.y[i];
            r1.z[k]=x.z[i];

            r2.x[k]=y.x[j];
            r2.y[k]=y.y[j];
            r2.z[k]=y.z[j];

            xtm.x[k]=x.x[i];
            xtm.y[k]=x.y[i];
            xtm.z[k]=x.z[i];

            ytm.x[k]=y.x[j];
            ytm.y[k]=y.y[j];
            ytm.z[k]=y.z[j];

            k++;
        }
        else if(i!=-1) {
            BasicFunction::PrintErrorAndQuit("Wrong map!\n");
        }
    }
    KabschFast(r1, r2, k, &rms, t, u, mem);

    //evaluate score
    double di;
    double d02=d0*d0;

    int n_ali=k;
    float xrot[3];
    tmscore=0;
    for(k=0; k<n_ali; k++)
    {
        BasicFunction::transform(t, u, xtm.x[k], xtm.y[k], xtm.z[k], xrot[0], xrot[1], xrot[2]);
        di=BasicFunction::dist(xrot[0], xrot[1], xrot[2],
                ytm.x[k], ytm.y[k], ytm.z[k]);
        tmscore += 1/(1+di/d02);
    }
//    if(tmscore/Lnorm  > 0.5){
//    cout << "di: " << di << '\n' << "tmscore[k]: " << tmscore/Lnorm << endl;
//    }
    return tmscore/Lnorm; // no need to normalize this score because it will not be used for latter scoring
}


//Comprehensive TMscore search engine
// input:   two vector sets: x, y
//          an alignment invmap0[] between x and y
//          simplify_step: 1 or 40 or other integers
//          score_sum_method: 0 for score over all pairs
//                            8 for socre over the pairs with BasicFunction::dist<score_d8
// output:  the best rotaion matrix t, u that results in highest TMscore
double detailed_search(Coordinates &r1, Coordinates &r2, Coordinates &xtm, Coordinates &ytm,
                       Coordinates &xt, const Coordinates &x, const Coordinates &y, int ylen,
                       int invmap0[], float t[3], float u[3][3], int simplify_step,
                       float local_d0_search, float Lnorm, float score_d8, float d0, float * mem)
{
    //x is model, y is template, try to superpose onto y
    int i, j, k;
    float tmscore;
    float rmsd;

    k=0;
    for(i=0; i<ylen; i++)
    {
        j=invmap0[i];
        if(j>=0) //aligned
        {
            xtm.x[k]=x.x[j];
            xtm.y[k]=x.y[j];
            xtm.z[k]=x.z[j];

            ytm.x[k]=y.x[i];
            ytm.y[k]=y.y[i];
            ytm.z[k]=y.z[i];
            k++;
        }
    }

    //detailed search 40-->1
    tmscore = TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u, simplify_step,
                              &rmsd, local_d0_search, Lnorm, score_d8, d0, mem);
    return tmscore;
}

double detailed_search_standard( Coordinates &r1, Coordinates &r2,
                                 Coordinates &xtm, Coordinates &ytm, Coordinates &xt,
                                 const Coordinates &x, const Coordinates &y,
                                 int ylen, int invmap0[], float t[3], float u[3][3],
                                 int simplify_step, int score_sum_method, double local_d0_search,
                                 const bool& bNormalize, float Lnorm, float score_d8, float d0, float * mem)
{
    //x is model, y is template, try to superpose onto y
    int i, j, k;
    float tmscore;
    float rmsd;

    k=0;
    for(i=0; i<ylen; i++)
    {
        j=invmap0[i];
        if(j>=0) //aligned
        {
            xtm.x[k]=x.x[j];
            xtm.y[k]=x.y[j];
            xtm.z[k]=x.z[j];

            ytm.x[k]=y.x[i];
            ytm.y[k]=y.y[i];
            ytm.z[k]=y.z[i];
            k++;
        }
    }

    //detailed search 40-->1
    tmscore = TMscore8_search_standard( r1, r2, xtm, ytm, xt, k, t, u,
                                        simplify_step, score_sum_method, &rmsd, local_d0_search, score_d8, d0, mem);
    if (bNormalize)// "-i", to use standard_TMscore, then bNormalize=true, else bNormalize=false;
        tmscore = tmscore * k / Lnorm;

    return tmscore;
}

//compute the score quickly in three iterations
double get_score_fast( Coordinates &r1, Coordinates &r2, Coordinates &xtm, Coordinates &ytm,
                       const Coordinates &x, const Coordinates &y, int ylen, int invmap[],
                       float d0, float d0_search, float t[3], float u[3][3], float * mem)
{
    float rms, tmscore, tmscore1, tmscore2;
    int i, j, k;

    k=0;
    for(j=0; j<ylen; j++)
    {
        i=invmap[j];
        if(i>=0)
        {
            r1.x[k]=x.x[i];
            r1.y[k]=x.y[i];
            r1.z[k]=x.z[i];

            r2.x[k]=y.x[j];
            r2.y[k]=y.y[j];
            r2.z[k]=y.z[j];

            xtm.x[k]=x.x[i];
            xtm.y[k]=x.y[i];
            xtm.z[k]=x.z[i];

            ytm.x[k]=y.x[j];
            ytm.y[k]=y.y[j];
            ytm.z[k]=y.z[j];

            k++;
        }
        else if(i!=-1) BasicFunction::PrintErrorAndQuit("Wrong map!\n");
    }
    KabschFast(r1, r2, k, &rms, t, u, mem);

    //evaluate score
    double di;
    const int len=k;
    double dis[len];
    double d00=d0_search;
    double d002=d00*d00;
    double d02=d0*d0;

    int n_ali=k;
    float xrot[3];
    tmscore=0;
    for(k=0; k<n_ali; k++)
    {
        BasicFunction::BasicFunction::transform(t, u, xtm.x[k], xtm.y[k], xtm.z[k], xrot[0], xrot[1], xrot[2]);
        di=BasicFunction::BasicFunction::dist(xrot[0], xrot[1], xrot[2],
                ytm.x[k], ytm.y[k], ytm.z[k]);
        dis[k]=di;
        tmscore += 1/(1+di/d02);
    }



    //second iteration
    float d002t=d002;
    while(1)
    {
        j=0;
        for(k=0; k<n_ali; k++)
        {
            if(dis[k]<=d002t)
            {
                r1.x[j]=xtm.x[k];
                r1.y[j]=xtm.y[k];
                r1.z[j]=xtm.z[k];

                r2.x[j]=ytm.x[k];
                r2.y[j]=ytm.y[k];
                r2.z[j]=ytm.z[k];

                j++;
            }
        }
        //there are not enough feasible pairs, relieve the threshold
        if(j<3 && n_ali>3) d002t += 0.5;
        else break;
    }

    if(n_ali!=j)
    {
        KabschFast(r1, r2, j, &rms, t, u, mem);
        tmscore1=0;
        for(k=0; k<n_ali; k++)
        {
            BasicFunction::BasicFunction::transform(t, u, xtm.x[k],  xtm.y[k], xtm.z[k], xrot[0], xrot[1], xrot[2]);
            di=BasicFunction::BasicFunction::dist(xrot[0], xrot[1], xrot[2],
                    ytm.x[k], ytm.y[k], ytm.z[k]);
            dis[k]=di;
            tmscore1 += 1/(1+di/d02);
        }

        //third iteration
        d002t=d002+1;

        while(1)
        {
            j=0;
            for(k=0; k<n_ali; k++)
            {
                if(dis[k]<=d002t)
                {
                    r1.x[j]=xtm.x[k];
                    r1.y[j]=xtm.y[k];
                    r1.z[j]=xtm.z[k];
                    r2.x[j]=ytm.x[k];
                    r2.y[j]=ytm.y[k];
                    r2.z[j]=ytm.z[k];
                    j++;
                }
            }
            //there are not enough feasible pairs, relieve the threshold
            if(j<3 && n_ali>3) d002t += 0.5;
            else break;
        }

        //evaluate the score
        KabschFast(r1, r2, j, &rms, t, u, mem);
        tmscore2=0;
        for(k=0; k<n_ali; k++)
        {
            BasicFunction::BasicFunction::transform(t, u, xtm.x[k],  xtm.y[k], xtm.z[k], xrot[0], xrot[1], xrot[2]);
            di=BasicFunction::BasicFunction::dist(xrot[0], xrot[1], xrot[2],
                    ytm.x[k], ytm.y[k], ytm.z[k]);
            tmscore2 += 1/(1+di/d02);
        }
    }
    else
    {
        tmscore1=tmscore;
        tmscore2=tmscore;
    }

    if(tmscore1>=tmscore) tmscore=tmscore1;
    if(tmscore2>=tmscore) tmscore=tmscore2;
    return tmscore; // no need to normalize this score because it will not be used for latter scoring
}


//perform gapless threading to find the best initial alignment
//input: x, y, xlen, ylen
//output: y2x0 stores the best alignment: e.g.,
//y2x0[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0
//the jth element in y is aligned to a gap in x if i==-1
double get_initial(Coordinates &r1, Coordinates &r2,
                   Coordinates &xtm, Coordinates &ytm,
                   const Coordinates &x, const Coordinates &y, int xlen, int ylen, int *y2x,
                   float d0, float d0_search, const bool fast_opt,
                   float t[3], float u[3][3], float * mem)
{
    int min_len=std::min(xlen, ylen);
    if(min_len<=5) BasicFunction::PrintErrorAndQuit("Sequence is too short <=5!\n");

    int min_ali= min_len/2;              //minimum size of considered fragment
    if(min_ali<=5)  min_ali=5;
    int n1, n2;
    n1 = -ylen+min_ali;
    n2 = xlen-min_ali;

    int i, j, k, k_best;
    float tmscore, tmscore_max=-1;

    k_best=n1;
    for(k=n1; k<=n2; k+=(fast_opt)?5:1)
    {
        //get the map
        for(j=0; j<ylen; j++)
        {
            i=j+k;
            if(i>=0 && i<xlen) y2x[j]=i;
            else y2x[j]=-1;
        }

        //evaluate the map quickly in three iterations
        //this is not real tmscore, it is used to evaluate the goodness of the initial alignment
        tmscore=get_score_fast(r1, r2, xtm, ytm,
                               x, y, ylen, y2x, d0,d0_search, t, u, mem);
        if(tmscore>=tmscore_max)
        {
            tmscore_max=tmscore;
            k_best=k;
        }
    }

    //extract the best map
    k=k_best;
    for(j=0; j<ylen; j++)
    {
        i=j+k;
        if(i>=0 && i<xlen) y2x[j]=i;
        else y2x[j]=-1;
    }

    return tmscore_max;
}

void smooth(int *sec, int len)
{
    int i, j;
    //smooth single  --x-- => -----
    for (i=2; i<len-2; i++)
    {
        if(sec[i]==2 || sec[i]==4)
        {
            j=sec[i];
            if (sec[i-2]!=j && sec[i-1]!=j && sec[i+1]!=j && sec[i+2]!=j)
                sec[i]=1;
        }
    }

    //   smooth double
    //   --xx-- => ------
    for (i=0; i<len-5; i++)
    {
        //helix
        if (sec[i]!=2   && sec[i+1]!=2 && sec[i+2]==2 && sec[i+3]==2 &&
            sec[i+4]!=2 && sec[i+5]!= 2)
        {
            sec[i+2]=1;
            sec[i+3]=1;
        }

        //beta
        if (sec[i]!=4   && sec[i+1]!=4 && sec[i+2]==4 && sec[i+3]==4 &&
            sec[i+4]!=4 && sec[i+5]!= 4)
        {
            sec[i+2]=1;
            sec[i+3]=1;
        }
    }

    //smooth connect
    for (i=0; i<len-2; i++)
    {
        if (sec[i]==2 && sec[i+1]!=2 && sec[i+2]==2) sec[i+1]=2;
        else if(sec[i]==4 && sec[i+1]!=4 && sec[i+2]==4) sec[i+1]=4;
    }

}

int sec_str(float dis13, float dis14, float dis15,
            float dis24, float dis25, float dis35)
{
    int s=1;

    float delta=2.1;
    if (fabsf(dis15-6.37)<delta && fabsf(dis14-5.18)<delta &&
        fabsf(dis25-5.18)<delta && fabsf(dis13-5.45)<delta &&
        fabsf(dis24-5.45)<delta && fabsf(dis35-5.45)<delta)
    {
        s=2; //helix
        return s;
    }

    delta=1.42;
    if (fabsf(dis15-13  )<delta && fabsf(dis14-10.4)<delta &&
        fabsf(dis25-10.4)<delta && fabsf(dis13-6.1 )<delta &&
        fabsf(dis24-6.1 )<delta && fabsf(dis35-6.1 )<delta)
    {
        s=4; //strand
        return s;
    }

    if (dis15 < 8) s=3; //turn
    return s;
}


//1->coil, 2->helix, 3->turn, 4->strand
void make_sec(Coordinates &x, int len, char *sec)
{
    int j1, j2, j3, j4, j5;
    float d13, d14, d15, d24, d25, d35;
    for(int i=0; i<len; i++)
    {
        sec[i]=1;
        j1=i-2;
        j2=i-1;
        j3=i;
        j4=i+1;
        j5=i+2;

        if(j1>=0 && j5<len)
        {
            d13=sqrt(BasicFunction::dist(x.x[j1], x.y[j1], x.z[j1], x.x[j3], x.y[j3], x.z[j3]));
            d14=sqrt(BasicFunction::dist(x.x[j1], x.y[j1], x.z[j1], x.x[j4], x.y[j4], x.z[j4]));
            d15=sqrt(BasicFunction::dist(x.x[j1], x.y[j1], x.z[j1], x.x[j5], x.y[j5], x.z[j5]));
            d24=sqrt(BasicFunction::dist(x.x[j2], x.y[j2], x.z[j2], x.x[j4], x.y[j4], x.z[j4]));
            d25=sqrt(BasicFunction::dist(x.x[j2], x.y[j2], x.z[j2], x.x[j5], x.y[j5], x.z[j5]));
            d35=sqrt(BasicFunction::dist(x.x[j3], x.y[j3], x.z[j3], x.x[j5], x.y[j5], x.z[j5]));
            sec[i]=sec_str(d13, d14, d15, d24, d25, d35);
        }
    }
}


//get initial alignment from secondary structure alignment
//input: x, y, xlen, ylen
//output: y2x stores the best alignment: e.g.,
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0
//the jth element in y is aligned to a gap in x if i==-1
void get_initial_ss(AffineNeedlemanWunsch *affineNW,
                    const char *secx, const char *secy, int xlen, int ylen, int *y2x)
{
    static const int ss_mat[] = {
            0, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1
    };

//    static const int Ori_CLESUM_WS[18*18]={
//     //      A   C   D    E       F       G      H    I       K      L   M   N   P   Q   R   S   T  V  X
//            73,  20,  13,  -17,  -25,  -20,  -6,  -45,  -31, -23, -19, -11,  -2,  10,  25,  35,  16,  0,
//            20,  51,   7,   13,   15,    7,  13,  -96,  -74, -57, -50, -12, -13, -11, -12,  42,  12,  0,
//            13,   7,  53,   21,    3,   20,  -4,  -77,  -56, -43, -33,   0, -12,  -5,   3,   4,  29,  0,
//            -17,  13,  21,   52,   22,   22, -31, -124, -105, -88, -81, -22, -49, -44, -42, -10,  14,  0,
//            -25,  15,   3,   22,   36,   26, -22, -127, -108, -93, -84, -21, -47, -43, -48,  -5,  -6,  0,
//            -20,   7,  20,   22,   26,   50,  -5, -107,  -88, -73, -69, -16, -33, -32, -30,   0,   3,  0,
//            -6,  13,  -4,  -31,  -22,   -5,  69,  -51,  -34, -21, -13,  29,  21,  -8,  -1,   5,   8,  0,
//            -45, -96, -77, -124, -127, -107, -51,   23,   18,  13,   5, -62,  -4, -34, -55, -60, -87,  0,
//            -31, -74, -56, -105, -108,  -88, -34,   18,   23,  16,  21, -41,   1, -11, -34, -49, -62,  0,
//            -23, -57, -43,  -88,  -93,  -73, -21,   13,   16,  37,  13, -32,  16,  -2, -24, -34, -44,  0,
//            -19, -50, -33,  -81,  -84,  -69, -13,    5,   21,  13,  49,  -1,  12,  28,   5, -36, -24,  0,
//            -11, -12,   0,  -22,  -21,  -16,  29,  -62,  -41, -32,  -1,  74,   5,   8,  -4, -12,  26,  0,
//            -2, -13, -12,  -49,  -47,  -33,  21,   -4,    1,  16,  12,   5,  61,   7,   5,   8,  -7,  0,
//            10, -11,  -5,  -44,  -43,  -32,  -8,  -34,  -11,  -2,  28,   8,   7,  90,  15,  -3,  32,  0,
//            25, -12,   3,  -42,  -48,  -30,  -1,  -55,  -34, -24,   5,  -4,   5,  15, 104,   4, -13,  0,
//            35,  42,   4,  -10,   -5,    0,   5,  -60,  -49, -34, -36, -12,   8,  -3,   4,  66,   7,  0,
//            16,  12,  29,   14,   -6,    3,   8,  -87,  -62, -44, -24,  26,  -7,  32, -13,   7,  90,  0,
//            0,   0,   0,    0,    0,    0,   0,    0,    0,   0,   0,   0,   0,   0,   0,   0,   0,  0};//R
// A   B   C    D    E    F   G    H    I   J   K   L   M   N   O   P   Q R

    static const int map[256] = {
            0, 1, 2, 3, 4
    };
    //      A   C   D    E       F       G      H    I       K      L   M   N   P   Q   R   S   T  V  X
//    static const int map[256] = {
//            23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//            23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//            23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//            23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//            23,  0, 20,  1,  2,  3, 4,  5,  6,  7, 23, 8, 9, 10,  11, 23,
//            12,  13,  14, 15, 16, 23, 17, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//            23,  0, 20,  1,  2,  3, 4,  5,  6,  7, 23, 8, 0, 10,  11, 23,
//            12,  13,  14, 15, 16, 23, 17, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//            23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//            23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//            23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//            23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//            23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//            23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//            23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//            23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//    };



    static const AffineNeedlemanWunsch::matrix_t matrix = {
            "ss",
            ss_mat,
            map,
            5,
            1,
            0
    };



//    NWDP_TM(score, path, val, secx, secy, xlen, ylen, -1.0, y2x);
//    for(size_t i = 0; i < ylen; i++){
//        if(y2x[i]!=-1)
//        std::cout << i << "\t" << y2x[i] << "\t" << (int) secy[i] << "\t" << (int)  secx[y2x[i]] << std::endl;
//    }
    std::fill(y2x, y2x+ylen, -1);
    AffineNeedlemanWunsch::profile_t *profile = affineNW->profile_create(secy, ylen, &matrix);
    affineNW->align(profile, ylen, (const unsigned char * ) secx, xlen,  100, 0, y2x);

//    for(size_t i = 0; i < ylen; i++){
//        if(y2x[i]!=-1)
//            std::cout << i << "\t" << y2x[i] << "\t" <<  secy[i] << "\t" <<   secx[y2x[i]] << std::endl;
//    }
//    std::cout << std::endl;
}


// get_initial5 in TMalign fortran, get_intial_local in TMalign c by yangji
//get initial alignment of local structure superposition
//input: x, y, xlen, ylen
//output: y2x stores the best alignment: e.g.,
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0
//the jth element in y is aligned to a gap in x if i==-1
bool get_initial5(AffineNeedlemanWunsch *affineNW,
                   Coordinates &r1, Coordinates &r2, Coordinates &xtm, Coordinates &ytm,
                   const Coordinates &x, const Coordinates &y, int xlen, int ylen, int *y2x,
                   float d0, float d0_search, const bool fast_opt, const float D0_MIN, float * mem)
{
    float GL, rmsd;
    float t[3];
    float u[3][3];

    float d01 = d0 + 1.5;
    if (d01 < D0_MIN) d01 = D0_MIN;
    float d02 = d01*d01;

    float GLmax = 0;
    int aL = std::min(xlen, ylen);
    int *invmap = new int[ylen + 1];

    // jump on sequence1-------------->
    int n_jump1 = 0;
    if (xlen > 250)
        n_jump1 = 45;
    else if (xlen > 200)
        n_jump1 = 35;
    else if (xlen > 150)
        n_jump1 = 25;
    else
        n_jump1 = 15;
    if (n_jump1 > (xlen / 3))
        n_jump1 = xlen / 3;

    // jump on sequence2-------------->
    int n_jump2 = 0;
    if (ylen > 250)
        n_jump2 = 45;
    else if (ylen > 200)
        n_jump2 = 35;
    else if (ylen > 150)
        n_jump2 = 25;
    else
        n_jump2 = 15;
    if (n_jump2 > (ylen / 3))
        n_jump2 = ylen / 3;

    // fragment to superimpose-------------->
    int n_frag[2] = { 20, 100 };
    if (n_frag[0] > (aL / 3))
        n_frag[0] = aL / 3;
    if (n_frag[1] > (aL / 2))
        n_frag[1] = aL / 2;

    // start superimpose search-------------->
    if (fast_opt)
    {
        n_jump1*=5;
        n_jump2*=5;
    }
    bool flag = false;
    for (int i_frag = 0; i_frag < 2; i_frag++)
    {
        int m1 = xlen - n_frag[i_frag] + 1;
        int m2 = ylen - n_frag[i_frag] + 1;

        for (int i = 0; i<m1; i = i + n_jump1) //index starts from 0, different from FORTRAN
        {
            for (int j = 0; j<m2; j = j + n_jump2)
            {
                for (int k = 0; k<n_frag[i_frag]; k++) //fragment in y
                {
                    r1.x[k] = x.x[k + i];
                    r1.y[k] = x.y[k + i];
                    r1.z[k] = x.z[k + i];

                    r2.x[k] = y.x[k + j];
                    r2.y[k] = y.y[k + j];
                    r2.z[k] = y.z[k + j];
                }

                // superpose the two structures and rotate it
                KabschFast(r1, r2, n_frag[i_frag], &rmsd, t, u, mem);

                float gap_open = 0.0;
                //NWDP_TM(score, path, val,
                //        x, y, xlen, ylen, t, u, d02, gap_open, invmap, mem);
                std::fill(invmap, invmap+ylen, -1);
                AffineNeedlemanWunsch::profile_t *profile = affineNW->profile_xyz_create(NULL, ylen, y.x, y.y, y.z);
                affineNW->alignXYZ(profile, ylen, xlen, x.x, x.y, x.z,
                                                                              d02, t, u, gap_open, 0.0, invmap);

                GL = get_score_fast(r1, r2, xtm, ytm, x, y, ylen,
                                    invmap, d0, d0_search, t, u, mem);
                if (GL>GLmax)
                {
                    GLmax = GL;
                    for (int ii = 0; ii<ylen; ii++) y2x[ii] = invmap[ii];
                    flag = true;
                }
            }
        }
    }

    delete[] invmap;
    return flag;
}

void score_matrix_rmsd_sec( Coordinates &r1,  Coordinates &r2,
                            float **score, const char *secx, const char *secy,
                            const Coordinates &x, const Coordinates &y, int xlen, int ylen,
                            int *y2x, const float D0_MIN, float d0, float * mem)
{
    float t[3], u[3][3];
    float rmsd, dij;
    float d01=d0+1.5;
    if(d01 < D0_MIN) d01=D0_MIN;
    float d02=d01*d01;

    float xx[3];
    int i, k=0;
    for(int j=0; j<ylen; j++)
    {
        i=y2x[j];
        if(i>=0)
        {
            r1.x[k]=x.x[i];
            r1.y[k]=x.y[i];
            r1.z[k]=x.z[i];

            r2.x[k]=y.x[j];
            r2.y[k]=y.y[j];
            r2.z[k]=y.z[j];

            k++;
        }
    }
    KabschFast(r1, r2, k, &rmsd, t, u, mem);


    for(int ii=0; ii<xlen; ii++)
    {
        BasicFunction::BasicFunction::transform(t, u, x.x[ii], x.y[ii], x.z[ii], xx[0], xx[1], xx[2]);
        for(int jj=0; jj<ylen; jj++)
        {
            dij=BasicFunction::BasicFunction::dist(xx[0], xx[1], xx[2], y.x[jj], y.y[jj], y.z[jj]);
            if (secx[ii]==secy[jj])
                score[ii+1][jj+1] = 1.0/(1+dij/d02) + 0.5;
            else
                score[ii+1][jj+1] = 1.0/(1+dij/d02);
        }
    }
}


//get initial alignment from secondary structure and previous alignments
//input: x, y, xlen, ylen
//output: y2x stores the best alignment: e.g.,
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0
//the jth element in y is aligned to a gap in x if i==-1
void get_initial_ssplus(Coordinates &r1, Coordinates &r2, float **score, bool **path,
                        float **val, const char *secx, const char *secy,
                        const Coordinates &x, const Coordinates &y,
                        int xlen, int ylen, int *y2x0, int *y2x, const double D0_MIN, double d0, float * mem)
{
    //create score matrix for DP
    score_matrix_rmsd_sec(r1, r2, score, secx, secy, x, y, xlen, ylen,
                          y2x0, D0_MIN,d0, mem);

    float gap_open=-1.0;
    NWDP_TM(score, path, val, xlen, ylen, gap_open, y2x);
}

void find_max_frag(const Coordinates &x, int len, int *start_max,
                   int *end_max, float dcu0, const bool fast_opt)
{
    int r_min, fra_min=4;           //minimum fragment for search
    if (fast_opt) fra_min=8;
    int start;
    int Lfr_max=0;

    r_min= (int) (len*1.0/3.0); //minimum fragment, in case too small protein
    if(r_min > fra_min) r_min=fra_min;

    int inc=0;
    float dcu0_cut=dcu0*dcu0;;
    float dcu_cut=dcu0_cut;

    while(Lfr_max < r_min)
    {
        Lfr_max=0;
        int j=1;    //number of residues at nf-fragment
        start=0;
        for(int i=1; i<len; i++)
        {
            if(BasicFunction::dist(x.x[i-1], x.y[i-1], x.z[i-1], x.x[i], x.y[i], x.z[i]) < dcu_cut)
            {
                j++;

                if(i==(len-1))
                {
                    if(j > Lfr_max)
                    {
                        Lfr_max=j;
                        *start_max=start;
                        *end_max=i;
                    }
                    j=1;
                }
            }
            else
            {
                if(j>Lfr_max)
                {
                    Lfr_max=j;
                    *start_max=start;
                    *end_max=i-1;
                }

                j=1;
                start=i;
            }
        }// for i;

        if(Lfr_max < r_min)
        {
            inc++;
            double dinc=pow(1.1, (double) inc) * dcu0;
            dcu_cut= dinc*dinc;
        }
    }//while <;
}

//perform fragment gapless threading to find the best initial alignment
//input: x, y, xlen, ylen
//output: y2x0 stores the best alignment: e.g.,
//y2x0[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0
//the jth element in y is aligned to a gap in x if i==-1
double get_initial_fgt(Coordinates &r1, Coordinates &r2, Coordinates &xtm, Coordinates &ytm,
                       const Coordinates &x, const Coordinates &y, int xlen, int ylen,
                       int *y2x, float d0, float d0_search,
                       float dcu0, const bool fast_opt, float t[3], float u[3][3],
                       float * mem)
{
    int fra_min=4;           //minimum fragment for search
    if (fast_opt) fra_min=8;
    int fra_min1=fra_min-1;  //cutoff for shift, save time

    int xstart=0, ystart=0, xend=0, yend=0;

    find_max_frag(x, xlen, &xstart, &xend, dcu0, fast_opt);
    find_max_frag(y, ylen, &ystart, &yend, dcu0, fast_opt);


    int Lx = xend-xstart+1;
    int Ly = yend-ystart+1;
    int *ifr, *y2x_;
    int L_fr=std::min(Lx, Ly);
    ifr= new int[L_fr];
    y2x_= new int[ylen+1];

    //select what piece will be used (this may araise ansysmetry, but
    //only when L1=L2 and Lfr1=Lfr2 and L1 ne Lfr1
    //if L1=Lfr1 and L2=Lfr2 (normal proteins), it will be the same as initial1

    if(Lx<Ly || (Lx==Ly && xlen<=ylen))
    {
        for(int i=0; i<L_fr; i++) ifr[i]=xstart+i;
    }
    else if(Lx>Ly || (Lx==Ly && xlen>ylen))
    {
        for(int i=0; i<L_fr; i++) ifr[i]=ystart+i;
    }


    int L0=std::min(xlen, ylen); //non-redundant to get_initial1
    if(L_fr==L0)
    {
        int n1= (int)(L0*0.1); //my index starts from 0
        int n2= (int)(L0*0.89);

        int j=0;
        for(int i=n1; i<= n2; i++)
        {
            ifr[j]=ifr[i];
            j++;
        }
        L_fr=j;
    }


    //gapless threading for the extracted fragment
    double tmscore, tmscore_max=-1;

    if(Lx<Ly || (Lx==Ly && xlen<=ylen))
    {
        int L1=L_fr;
        int min_len=std::min(L1, ylen);
        int min_ali= (int) (min_len/2.5);              //minimum size of considered fragment
        if(min_ali<=fra_min1)  min_ali=fra_min1;
        int n1, n2;
        n1 = -ylen+min_ali;
        n2 = L1-min_ali;

        int i, j, k;
        for(k=n1; k<=n2; k+=(fast_opt)?3:1)
        {
            //get the map
            for(j=0; j<ylen; j++)
            {
                i=j+k;
                if(i>=0 && i<L1) y2x_[j]=ifr[i];
                else             y2x_[j]=-1;
            }

            //evaluate the map quickly in three iterations
            tmscore=get_score_fast(r1, r2, xtm, ytm, x, y, ylen, y2x_,
                                   d0, d0_search, t, u, mem);

            if(tmscore>=tmscore_max)
            {
                tmscore_max=tmscore;
                for(j=0; j<ylen; j++) y2x[j]=y2x_[j];
            }
        }
    }
    else
    {
        int L2=L_fr;
        int min_len=std::min(xlen, L2);
        int min_ali= (int) (min_len/2.5);              //minimum size of considered fragment
        if(min_ali<=fra_min1)  min_ali=fra_min1;
        int n1, n2;
        n1 = -L2+min_ali;
        n2 = xlen-min_ali;

        int i, j, k;

        for(k=n1; k<=n2; k++)
        {
            //get the map
            for(j=0; j<ylen; j++) y2x_[j]=-1;

            for(j=0; j<L2; j++)
            {
                i=j+k;
                if(i>=0 && i<xlen) y2x_[ifr[j]]=i;
            }

            //evaluate the map quickly in three iterations
            tmscore=get_score_fast(r1, r2, xtm, ytm,
                                   x, y, ylen, y2x_, d0,d0_search, t, u, mem);
            if(tmscore>=tmscore_max)
            {
                tmscore_max=tmscore;
                for(j=0; j<ylen; j++) y2x[j]=y2x_[j];
            }
        }
    }


    delete [] ifr;
    delete [] y2x_;
    return tmscore_max;
}

//heuristic run of dynamic programing iteratively to find the best alignment
//input: initial rotation matrix t, u
//       vectors x and y, d0
//output: best alignment that maximizes the TMscore, will be stored in invmap
double DP_iter(AffineNeedlemanWunsch * affineNW,
               Coordinates &r1, Coordinates &r2,
               Coordinates &xtm, Coordinates &ytm,
               Coordinates &xt, const Coordinates &x, const Coordinates &y, int xlen, int ylen, float t[3], float u[3][3],
               int invmap0[], int g1, int g2, int iteration_max, float local_d0_search,
               float Lnorm, float d0, float score_d8, float * mem)
{
    float gap_open[2]={-0.6, 0};
    float rmsd;
    int *invmap=new int[ylen+1];

    int iteration, i, j, k;
    float tmscore, tmscore_max, tmscore_old=0;
    int simplify_step=40;
    tmscore_max=-1;

    //double d01=d0+1.5;
    float d02=d0*d0;
    for(int g=g1; g<g2; g++)
    {
        for(iteration=0; iteration<iteration_max; iteration++)
        {
//            NWDP_TM(score, path, val, x, y, xlen, ylen,
//                    t, u, d02, gap_open[g], invmap, mem);
            std::fill(invmap, invmap+ylen, -1);
            AffineNeedlemanWunsch::profile_t *profile = affineNW->profile_xyz_create(NULL, ylen, y.x, y.y, y.z);
            affineNW->alignXYZ(profile, ylen, xlen, x.x, x.y, x.z,
                                                                          d02, t, u, -gap_open[g], 0.0, invmap);
//            std::cout << result.start_query << "\t" << result.start_target << std::endl;
//            std::cout << result.end_query << "\t" << result.end_target << std::endl;


            k=0;
            for(j=0; j<ylen; j++)
            {
                i=invmap[j];
                if(i>=0) //aligned
                {
                    xtm.x[k]=x.x[i];
                    xtm.y[k]=x.y[i];
                    xtm.z[k]=x.z[i];

                    ytm.x[k]=y.x[j];
                    ytm.y[k]=y.y[j];
                    ytm.z[k]=y.z[j];
                    k++;
                }
            }

            tmscore = TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u,
                                      simplify_step, &rmsd, local_d0_search,
                                      Lnorm, score_d8, d0, mem);


            if(tmscore>tmscore_max)
            {
                tmscore_max=tmscore;
                for(i=0; i<ylen; i++) invmap0[i]=invmap[i];
            }

            if(iteration>0)
            {
                if(fabs(tmscore_old-tmscore)<0.000001) break;
            }
            tmscore_old=tmscore;
        }// for iteration

    }//for gapopen


    delete []invmap;
    return tmscore_max;
}


double standard_TMscore(Coordinates &r1, Coordinates &r2, Coordinates &xtm, Coordinates &ytm,
                        Coordinates &xt, Coordinates &x, Coordinates &y, int ylen, int invmap[],
                        int& L_ali, float& RMSD, float D0_MIN, float Lnorm, float d0,
                        float score_d8, float t[3], float u[3][3], float * mem)
{
    D0_MIN = 0.5;
    Lnorm = ylen;
    if (Lnorm > 21)
        d0 = (1.24*pow((Lnorm*1.0 - 15), 1.0 / 3) - 1.8);
    else
        d0 = D0_MIN;
    if (d0 < D0_MIN)
        d0 = D0_MIN;
    double d0_input = d0;// Scaled by seq_min

    double tmscore;// collected alined residues from invmap
    int n_al = 0;
    int i;
    for (int j = 0; j<ylen; j++)
    {
        i = invmap[j];
        if (i >= 0)
        {
            xtm.x[n_al] = x.x[i];
            xtm.y[n_al] = x.y[i];
            xtm.z[n_al] = x.z[i];
//            std::cout << "x: " << i << " " << j << " " << x.x[i] << " " << x.y[i] << " " << x.z[i] << std::endl;
            ytm.x[n_al] = y.x[j];
            ytm.y[n_al] = y.y[j];
            ytm.z[n_al] = y.z[j];
//            std::cout << "y: " << i << " " << j << " " << y.x[j] << " " << y.y[j] << " " << y.z[j] << std::endl;

            r1.x[n_al] = x.x[i];
            r1.y[n_al] = x.y[i];
            r1.z[n_al] = x.z[i];

            r2.x[n_al] = y.x[j];
            r2.y[n_al] = y.y[j];
            r2.z[n_al] = y.z[j];
            n_al++;
        }
        else if (i != -1) {
            BasicFunction::PrintErrorAndQuit("Wrong map!\n");
        }
    }
    L_ali = n_al;

    KabschFast(r1, r2, n_al, &RMSD, t, u, mem);
    if(u[0][0]==NAN){
        return 0;
    }
    RMSD = sqrt( RMSD/(1.0*n_al) );

    int temp_simplify_step = 40;
    int temp_score_sum_method = 0;
    float rms = 0.0;
    tmscore = TMscore8_search_standard(r1, r2, xtm, ytm, xt, n_al, t, u,
                                       temp_simplify_step, temp_score_sum_method, &rms, d0_input,
                                       score_d8, d0, mem);
    tmscore = tmscore * n_al / (1.0*Lnorm);

    return tmscore;
}

template <class A> void DeleteArray(A *** array, int Narray)
{
    for(int i=0; i<Narray; i++)
        if(*(*array+i)) delete [] *(*array+i);
    if(Narray) delete [] (*array);
    (*array)=NULL;
}

template <class A> void NewArray(A *** array, int Narray1, int Narray2)
{
    *array=new A* [Narray1];
    for(int i=0; i<Narray1; i++) *(*array+i)=new A [Narray2];
}

/* entry function for TMalign */
int TMalign_main(
        AffineNeedlemanWunsch * affineNW,
        const Coordinates &xa, const Coordinates &ya,
        const char *seqx, const char *seqy, const char *secx, const char *secy,
        float t0[3], float u0[3][3],
        float &TM1, float &TM2, float &TM3, float &TM4, float &TM5,
        float &d0_0, float &TM_0,
        float &d0A, float &d0B, float &d0u, float &d0a, float &d0_out,
        string &seqM, string &seqxA, string &seqyA,
        float &rmsd0, float &Liden, int &n_ali, int &n_ali8,
        const int xlen, const int ylen, const float Lnorm_ass,
        const float d0_scale, const bool I_opt, const bool a_opt,
        const bool u_opt, const bool d_opt, const bool fast_opt, float * mem)
{
    int minlen = min(xlen, ylen);

    float D0_MIN;        //for d0
    float Lnorm;         //normalization length
    float score_d8,d0,d0_search,dcu0;//for TMscore search
    float t[3], u[3][3]; //Kabsch translation vector and rotation matrix
    float **score;     // Input score table for dynamic programming
    bool   **path;     // for dynamic programming
    float **val;       // for dynamic programming
    Coordinates xtm(minlen);
    Coordinates ytm(minlen);
    Coordinates xt(xlen);
    Coordinates r1(minlen);
    Coordinates r2(minlen);


    /***********************/
    /* allocate memory     */
    /***********************/
    NewArray(&score, xlen+1, ylen+1);
    NewArray(&path, xlen+1, ylen+1);
    NewArray(&val, xlen+1, ylen+1);

//    NewArray(&xtm, minlen, 3);
//    NewArray(&ytm, minlen, 3);
//    NewArray(&xt, xlen, 3);
//    NewArray(&r1, minlen, 3);
//    NewArray(&r2, minlen, 3);

    /***********************/
    /*    parameter set    */
    /***********************/
    parameter_set4search(xlen, ylen, D0_MIN, Lnorm,
                         score_d8, d0, d0_search, dcu0);
    int simplify_step    = 40; //for similified search engine
    int score_sum_method = 8;  //for scoring method, whether only sum over pairs with dis<score_d8

    int i;
    int *invmap0         = new int[ylen+1];
    int *invmap          = new int[ylen+1];
    float TM, TMmax=-1;
    for(i=0; i<ylen; i++) invmap0[i]=-1;

    float ddcc=0.4;
    if (Lnorm <= 40) ddcc=0.1;   //Lnorm was setted in parameter_set4search
    float local_d0_search = d0_search;


    /******************************************************/
    /*    get initial alignment with gapless threading    */
    /******************************************************/

    get_initial(r1, r2, xtm, ytm, xa, ya, xlen, ylen, invmap0, d0,
                d0_search, fast_opt, t, u, mem);
    TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, ylen, invmap0,
                         t, u, simplify_step, local_d0_search, Lnorm,
                         score_d8, d0, mem);
    if (TM>TMmax) TMmax = TM;
    //run dynamic programing iteratively to find the best alignment
    TM = DP_iter(affineNW, r1, r2, xtm, ytm, xt, xa, ya,
                  xlen, ylen, t, u, invmap, 0, 2, (fast_opt)?2:30, local_d0_search,
                  Lnorm, d0, score_d8, mem);
    if (TM>TMmax)
    {
        TMmax = TM;
        for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
    }


    /************************************************************/
    /*    get initial alignment based on secondary structure    */
    /************************************************************/
    get_initial_ss(affineNW, secx, secy, xlen, ylen, invmap);
    TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, ylen, invmap,
                         t, u, simplify_step, local_d0_search, Lnorm,
                         score_d8, d0, mem);
    if (TM>TMmax)
    {
        TMmax = TM;
        for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
    }
    if (TM > TMmax*0.2)
    {
        TM = DP_iter(affineNW, r1, r2, xtm, ytm, xt, xa, ya,
                     xlen, ylen, t, u, invmap, 0, 2, (fast_opt)?2:30,
                     local_d0_search, Lnorm, d0, score_d8, mem);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
    }


    /************************************************************/
    /*    get initial alignment based on local superposition    */
    /************************************************************/
    //=initial5 in original TM-align
    if (get_initial5(affineNW, r1, r2, xtm, ytm, xa, ya,
                      xlen, ylen, invmap, d0, d0_search, fast_opt, D0_MIN, mem))
    {
        TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, ylen,
                             invmap, t, u, simplify_step,
                             local_d0_search, Lnorm, score_d8, d0, mem);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
        if (TM > TMmax*ddcc)
        {
            TM = DP_iter(affineNW, r1, r2, xtm, ytm, xt, xa, ya,
                         xlen, ylen, t, u, invmap, 0, 2, 2, local_d0_search,
                         Lnorm, d0, score_d8, mem);
            if (TM>TMmax)
            {
                TMmax = TM;
                for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
            }
        }
    }
    else
        cerr << "\n\nWarning: initial alignment from local superposition fail!\n\n" << endl;


    /********************************************************************/
    /* get initial alignment by local superposition+secondary structure */
    /********************************************************************/
    //=initial3 in original TM-align
    get_initial_ssplus(r1, r2, score, path, val, secx, secy, xa, ya,
                       xlen, ylen, invmap0, invmap, D0_MIN, d0, mem);
    TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, ylen, invmap,
                         t, u, simplify_step,  local_d0_search, Lnorm,
                         score_d8, d0, mem);
    if (TM>TMmax)
    {
        TMmax = TM;
        for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
    }
    if (TM > TMmax*ddcc)
    {
        TM = DP_iter(affineNW, r1, r2, xtm, ytm, xt, xa, ya,
                     xlen, ylen, t, u, invmap, 0, 2, (fast_opt)?2:30,
                     local_d0_search, Lnorm, d0, score_d8, mem);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
    }


    /*******************************************************************/
    /*    get initial alignment based on fragment gapless threading    */
    /*******************************************************************/
    //=initial4 in original TM-align
    //TODO
    get_initial_fgt(r1, r2, xtm, ytm, xa, ya, xlen, ylen,
                    invmap, d0, d0_search, dcu0, fast_opt, t, u, mem);
    TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, ylen, invmap,
                         t, u, simplify_step, local_d0_search, Lnorm,
                         score_d8, d0, mem);
    if (TM>TMmax)
    {
        TMmax = TM;
        for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
    }
    if (TM > TMmax*ddcc)
    {
        TM = DP_iter(affineNW, r1, r2, xtm, ytm, xt, xa, ya,
                     xlen, ylen, t, u, invmap, 1, 2, 2, local_d0_search,
                     Lnorm, d0, score_d8, mem);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
    }

    //*******************************************************************//
    //    The alignment will not be changed any more in the following    //
    //*******************************************************************//
    //check if the initial alignment is generated approately
    bool flag=false;
    for(i=0; i<ylen; i++)
    {
        if(invmap0[i]>=0)
        {
            flag=true;
            break;
        }
    }
    if(!flag)
    {
        cout << "There is no alignment between the two proteins!" << endl;
        cout << "Program stop with no result!" << endl;
        return 1;
    }


    //********************************************************************//
    //    Detailed TMscore search engine --> prepare for final TMscore    //
    //********************************************************************//
    //run detailed TMscore search engine for the best alignment, and
    //extract the best rotation matrix (t, u) for the best alginment
    simplify_step=1;
    if (fast_opt) simplify_step=40;
    score_sum_method=8;
    TM = detailed_search_standard(r1, r2, xtm, ytm, xt, xa, ya, ylen,
                                  invmap0, t, u, simplify_step, score_sum_method, local_d0_search,
                                  false, Lnorm, score_d8, d0, mem);

    //select pairs with dis<d8 for final TMscore computation and output alignment
    int k=0;
    int *m1, *m2;
    double d;
    m1=new int[xlen]; //alignd index in x
    m2=new int[ylen]; //alignd index in y
    BasicFunction::do_rotation(xa, xt, xlen, t, u);
    k=0;
    for(int j=0; j<ylen; j++)
    {
        i=invmap0[j];
        if(i>=0)//aligned
        {
            n_ali++;
            d =sqrt( BasicFunction::dist(xt.x[i], xt.y[i], xt.z[i], ya.x[j], ya.y[j], ya.z[j]));

            if (d <= score_d8 || (I_opt == true))
            {
                m1[k]=i;
                m2[k]=j;

                xtm.x[k]=xa.x[i];
                xtm.y[k]=xa.y[i];
                xtm.z[k]=xa.z[i];
                ytm.x[k]=ya.x[j];
                ytm.y[k]=ya.y[j];
                ytm.z[k]=ya.z[j];
                r1.x[k] = xt.x[i];
                r1.y[k] = xt.y[i];
                r1.z[k] = xt.z[i];
                r2.x[k] = ya.x[j];
                r2.y[k] = ya.y[j];
                r2.z[k] = ya.z[j];
                k++;
            }
        }
    }
    n_ali8=k;

    KabschFast(r1, r2, n_ali8, &rmsd0, t, u, mem);// rmsd0 is used for final output, only recalculate rmsd0, not t & u
    rmsd0 = sqrt(rmsd0 / n_ali8);


    //****************************************//
    //              Final TMscore             //
    //    Please set parameters for output    //
    //****************************************//
    float rmsd;
    simplify_step=1;
    score_sum_method=0;
    float Lnorm_0=ylen;


    //normalized by length of structure A
    parameter_set4final(Lnorm_0, D0_MIN, Lnorm,
                        d0, d0_search);
    d0A=d0;
    d0_0=d0A;
    local_d0_search = d0_search;
    TM1 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0, simplify_step,
                          &rmsd, local_d0_search, Lnorm, score_d8, d0, mem);
    TM_0 = TM1;

    //normalized by length of structure B
    parameter_set4final(xlen+0.0, D0_MIN, Lnorm,
                        d0, d0_search);
    d0B=d0;
    local_d0_search = d0_search;
    TM2 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t, u, simplify_step,
                          &rmsd, local_d0_search, Lnorm, score_d8, d0, mem);

    if (a_opt)
    {
        //normalized by average length of structures A, B
        Lnorm_0=(xlen+ylen)*0.5;
        parameter_set4final(Lnorm_0, D0_MIN, Lnorm,
                            d0, d0_search);
        d0a=d0;
        d0_0=d0a;
        local_d0_search = d0_search;

        TM3 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
                              simplify_step, &rmsd, local_d0_search, Lnorm,
                              score_d8, d0, mem);
        TM_0=TM3;
    }
    if (u_opt)
    {
        //normalized by user assigned length
        parameter_set4final(Lnorm_ass, D0_MIN, Lnorm,
                            d0, d0_search);
        d0u=d0;
        d0_0=d0u;
        Lnorm_0=Lnorm_ass;
        local_d0_search = d0_search;
        TM4 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
                              simplify_step, &rmsd, local_d0_search, Lnorm,
                              score_d8, d0, mem);
        TM_0=TM4;
    }
    if (d_opt)
    {
        //scaled by user assigned d0
        parameter_set4scale(ylen, d0_scale, Lnorm,
                            d0, d0_search);
        d0_out=d0_scale;
        d0_0=d0_scale;
        //Lnorm_0=ylen;
        local_d0_search = d0_search;
        TM5 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
                              simplify_step, &rmsd, local_d0_search, Lnorm,
                              score_d8, d0, mem);
        TM_0=TM5;
    }

    /* derive alignment from superposition */
    int ali_len=xlen+ylen; //maximum length of alignment
    seqxA.assign(ali_len,'-');
    seqM.assign( ali_len,' ');
    seqyA.assign(ali_len,'-');

    BasicFunction::do_rotation(xa, xt, xlen, t, u);

    int kk=0, i_old=0, j_old=0;
    d=0;
    for(int k=0; k<n_ali8; k++)
    {
        for(int i=i_old; i<m1[k]; i++)
        {
            //align x to gap
            seqxA[kk]=seqx[i];
            seqyA[kk]='-';
            seqM[kk]=' ';
            kk++;
        }

        for(int j=j_old; j<m2[k]; j++)
        {
            //align y to gap
            seqxA[kk]='-';
            seqyA[kk]=seqy[j];
            seqM[kk]=' ';
            kk++;
        }

        seqxA[kk]=seqx[m1[k]];
        seqyA[kk]=seqy[m2[k]];
        Liden+=(seqxA[kk]==seqyA[kk]);
        d=sqrt(BasicFunction::dist(xt.x[m1[k]], xt.y[m1[k]], xt.z[m1[k]], ya.x[m2[k]], ya.y[m2[k]], ya.z[m2[k]]));

        if(d<d0_out) seqM[kk]=':';
        else         seqM[kk]='.';
        kk++;
        i_old=m1[k]+1;
        j_old=m2[k]+1;
    }

    //tail
    for(int i=i_old; i<xlen; i++)
    {
        //align x to gap
        seqxA[kk]=seqx[i];
        seqyA[kk]='-';
        seqM[kk]=' ';
        kk++;
    }
    for(int j=j_old; j<ylen; j++)
    {
        //align y to gap
        seqxA[kk]='-';
        seqyA[kk]=seqy[j];
        seqM[kk]=' ';
        kk++;
    }
    seqxA=seqxA.substr(0,kk);
    seqyA=seqyA.substr(0,kk);
    seqM =seqM.substr(0,kk);

    /* free memory */
    delete [] invmap0;
    delete [] invmap;
    delete [] m1;
    delete [] m2;
    DeleteArray(&score, xlen+1);
    DeleteArray(&path, xlen+1);
    DeleteArray(&val, xlen+1);

    free(xtm.x);
    free(xtm.y);
    free(xtm.z);

    free(ytm.x);
    free(ytm.y);
    free(ytm.z);

    free(xt.x);
    free(xt.y);
    free(xt.z);

    free(r1.x);
    free(r1.y);
    free(r1.z);

    free(r2.x);
    free(r2.y);
    free(r2.z);


//    DeleteArray(&xtm, minlen);
//    DeleteArray(&ytm, minlen);
//    DeleteArray(&xt, xlen);
//    DeleteArray(&r1, minlen);
//    DeleteArray(&r2, minlen);
    return 0; // zero for no exception
}
