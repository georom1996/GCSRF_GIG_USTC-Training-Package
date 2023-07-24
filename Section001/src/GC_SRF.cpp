//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//
// g++ -c SX_deconv.cpp && g++ SX_deconv.o sac_lpz.o -o SX_deconv
// rmean rtr taper
// bp c 0.033 0.125 p 2
// decimate 5
// decimate 2
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
typedef float SAMPLE; /* data type of seismogram samples */
typedef float REAL;   /* general floating point type */

#include "stdio.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdarg.h>
#include <stdlib.h>
#include <string>
#include <time.h>

#include "./sac_lpz.h"
#define SHC_EPSILON 1.0e-25
#define SHC_PI 3.14159265358979323846
#define SHC_RAD_TO_DEG (180.0 / SHC_PI)
#define TRACE void
#define BC_LINELTH 132
#define SHE_OFFSET 1700
#define SHE_SPECERROR 1800         /* special errors 1801..1899 */
#define SHE_ZWDW (SHE_OFFSET + 20) /* zero time window set */
#define SHC_EPSILON 1.0e-25
#define SQRT3_2 0.8660254
/* half of square root of 3 */
#define SQRT3 1.7320508
/* square root of 3 */
#define SQRT1_3 5.77350269e-1
#define SQRT1_2 7.07106781e-1
#define SQRT2_3 8.16496581e-1
//#define Abs(x) ((x)<0?-(x):(x))
#define Abs(x) ((x) < 0 ? -(x) : (x))
#define Nint(x) ((x) > 0 ? (int)((x) + 0.5) : (int)((x)-0.5))
extern "C" {
#include "./sac_lpz.h"
}

using namespace std;

// of the functions

void fold(float z[], float r[], long lth, float *zz, float *rr, float *zr);
void cut_data(float in_data[], float delta, float begin_time, float end_time,
              float out_data[]);
void get_spikefilter(float *start, int trclth, float ff[], float delta,
                     float reg);
void mt_levinson(float *r, float *g, float *f, int m);
void spiking(float *start, int trclth, float f[]);
void mt_fold(long la, float a[], long lb, float b[], float c[]);
void mt_rot2(SAMPLE *n, SAMPLE *e, long lth, REAL angle, SAMPLE r[],
             SAMPLE t[]);
float sc_polar2(SAMPLE zarr_ori[], SAMPLE rarr_ori[], float delta,
                float mdir_t1, float mdir_t2);
void sc_trccorr(REAL z[], REAL r[], long lth, REAL *zz, REAL *rr, REAL *zr);
void mt_rot3(SAMPLE *z, SAMPLE *n, SAMPLE *e, long lth, REAL azim, REAL inci,
             int type, SAMPLE l[], SAMPLE q[], SAMPLE t[]);
char *myFormatStringByFun(char *format, ...);

int main(int argc, char *argv[]) {
  clock_t start, end; //定义clock_t变量
  start = clock();    //开始时间
  string option_str;
  string name_Z = "";
  string name_N = "";
  string name_E = "";
  string spike_t1_str = "";
  string spike_t2_str = "";
  string reg_str = "";
  string inc_str = "";
  string detial_flag = "";
  float spike_t1 = 0.0;
  float spike_t2 = 0.0;
  float regnum = 3.0;
  float incnum = 0.0;

  for (int i = 0; i < argc; i++) {
    option_str = argv[i];
    if (option_str == "-Z") {
      name_Z = argv[i + 1];
      i++;
      continue;
    }
    if (option_str == "-N") {
      name_N = argv[i + 1];
      i++;
      continue;
    }
    if (option_str == "-E") {
      name_E = argv[i + 1];
      i++;
      continue;
    }
    if (option_str == "-t1") {

      spike_t1_str = argv[i + 1];
      spike_t1 = atof(spike_t1_str.c_str());
      i++;
      continue;
    }

    if (option_str == "-reg") {
      reg_str = argv[i + 1];
      regnum = atof(reg_str.c_str());
      i++;
      continue;
    }

    if (option_str == "-D") {
      // detial_flag=argv[i+1];
      detial_flag = "detail_info";
      i++;
      continue;
    }
  }

  // cout<<name_Z<<"||"<<name_R<<"||"<<detial_flag<<"||"<<endl;
  if (name_Z == "" && name_N == "" && detial_flag == "") {
    cout << "------------------------------------------------------------------"
            "----"
         << endl;
    cout << "+            Program for calculate deconvolution of GC_SRF        "
            "   +"
         << endl;
    cout << "+             Version 1.02 Date 2021/07/06 By GeoRom/zhou         "
            "   +"
         << endl;
    cout << "+                 Contact with :zhagnzhou3@gig.ac.cn              "
            "   +"
         << endl;
    cout << "------------------------------------------------------------------"
            "----"
         << endl;
    cout << "GC_SRF ERROR : Pleas check the options:" << endl;
    cout << "Right usage :" << endl;
    cout << "       GC_SRF -Z [file name of Z] -N [file name of N] -E [file "
            "name of E]  [-t1 t_spike_begin] -reg [control fil] -D [with "
            "internal files]"
         << endl;
    cout << "       -t1 for the begin of the spike" << endl;
    cout << "       -D for output the filter file" << endl;
    cout << "Output files: rf-[inci]-[win_len].[com].sac" << endl;
    return 1;
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // READING THE Z COMPONENT++++++++++++++++++++++++++++++
  float *data_Z;
  float delta_Z;
  long int npts_Z;
  SACHEAD header_Z;
  data_Z = read_sac(name_Z.c_str(), &header_Z);
  if (data_Z == 0) {
    cout << "!error when read Z sac data!" << endl;
    return -1;
  }
  delta_Z = header_Z.delta;
  npts_Z = header_Z.npts;
  // END READING THE Z COMPONENT++++++++++++++++++++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // READING THE N COMPONENT++++++++++++++++++++++++++++++
  float *data_N;
  float delta_N;
  long int npts_N;
  SACHEAD header_N;
  data_N = read_sac(name_N.c_str(), &header_N);
  if (data_N == 0) {
    cout << "!error when read N sac data!" << endl;
    return -1;
  }
  delta_N = header_N.delta;
  npts_N = header_N.npts;
  // cout << data_R[10000] << " " << delta_R << " " << npts_R << endl;
  // END READING THE N COMPONENT++++++++++++++++++++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // READING THE E COMPONENT++++++++++++++++++++++++++++++
  float *data_E;
  float delta_E;
  long int npts_E;
  SACHEAD header_E;
  data_E = read_sac(name_E.c_str(), &header_E);
  if (data_E == 0) {
    cout << "!error when read E sac data!" << endl;
    return -1;
  }
  delta_E = header_E.delta;
  npts_E = header_E.npts;
  // create L and Q components
  SACHEAD header_L = header_Z, header_Q = header_Z, header_T = header_Z;
  float delta_L = header_L.delta, delta_Q = header_Q.delta,
        delta_T = header_T.delta;
  long npts_L = header_L.npts, npts_Q = header_Q.npts, npts_T = header_T.npts;
  strcpy(header_L.kcmpnm, "L-COM");
  strcpy(header_Q.kcmpnm, "Q-COM");
  strcpy(header_T.kcmpnm, "T-COM");
  float azimnum = header_Z.baz;
  float data_L[int(npts_L)], data_Q[int(npts_Q)], data_T[int(npts_T)];
  // end of creating L and Q components

  // LOOP INCI AND LEN START
  for (incnum = 0; incnum <= 60; incnum = incnum + 4) {
    for (float win_len = 5; win_len <= 95; win_len = win_len + 5) {
      spike_t2 = spike_t1 + win_len;
      mt_rot3(data_Z, data_N, data_E, npts_Z, azimnum, incnum, 1, data_L,
              data_Q, data_T);
      // cout<<"ZNE2LQT PROCESSING: BAZ="<<azimnum<<" INCI="<<incnum<<"
      // WIN_LEN="<<win_len<<"    OK!"<<endl;
      cout << "Now for: INCI=" << incnum << " WIN_LEN=" << win_len << " OK!"
           << endl;
      // mt_rot2( data_Z, data_R, npts_Z,mdir_angle_zr2lq, data_L, data_Q );

      if (detial_flag != "") {
        write_sac("l-com.sac", header_L, data_L);
        write_sac("q-com.sac", header_Q, data_Q);
        write_sac("t-com.sac", header_T, data_T);
      }
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SACHEAD header_fil = header_Z, header_DecL = header_L,
              header_DecQ = header_Q, header_DecT = header_T;
      float spike_b = 190;
      float spike_e = 240;

      if (spike_t1 != 0.0) {
        spike_b = spike_t1;
      }
      spike_e = spike_t2;

      int indexlo = 0;
      int indexhi = 0;
      indexlo = int(spike_b / header_Q.delta) + 1;
      indexhi = int(spike_e / header_Q.delta) + 1;
      int mv_points = int(indexhi - indexlo + 1);
      int npts_fil = int(mv_points / 2);
      float filter_spk[npts_fil];
      float *mov;
      mov = data_Q;
      get_spikefilter(mov + indexlo, mv_points, filter_spk, header_Q.delta,
                      regnum);

      // cout<<spike_b<<"  "<<spike_e<<endl;

      header_fil.npts = npts_fil;
      header_fil.b = 0;
      strcpy(header_fil.kcmpnm, "SPK_FIL");
      if (detial_flag != "") {
        write_sac("fil.sac", header_fil, filter_spk);
      }

      float data_DecL[npts_L + npts_fil];
      float data_DecQ[npts_Q + npts_fil];
      float data_DecT[npts_T + npts_fil];
      float data_DecLN[npts_L + npts_fil];
      float data_DecQN[npts_Q + npts_fil];
      float data_DecTN[npts_T + npts_fil];
      mt_fold(npts_L, data_L, header_fil.npts, filter_spk, data_DecL);
      mt_fold(npts_Q, data_Q, header_fil.npts, filter_spk, data_DecQ);
      mt_fold(npts_T, data_T, header_fil.npts, filter_spk, data_DecT);

      int maxN = 0;
      float max = 0.0;
      for (int i = indexlo; i < indexhi; i++) {
        if (max < data_DecQ[i]) {
          max = data_DecQ[i];
          maxN = i;
        }
      }
      strcpy(header_DecL.kcmpnm, "DecL-OUT");
      strcpy(header_DecQ.kcmpnm, "DecQ-OUT");
      strcpy(header_DecT.kcmpnm, "DecT-OUT");
      header_DecL.npts = npts_L + header_fil.npts;
      header_DecQ.npts = npts_Q + header_fil.npts;
      header_DecT.npts = npts_T + header_fil.npts;

      for (int i = 0; i < header_DecL.npts; i++) {
        data_DecLN[i] = data_DecL[i] / max;
        data_DecQN[i] = data_DecQ[i] / max;
        data_DecTN[i] = data_DecT[i] / max;
      }
      header_DecL.b = -maxN * header_DecL.delta;
      header_DecQ.b = -maxN * header_DecQ.delta;
      header_DecT.b = -maxN * header_DecT.delta;
      header_DecL.t1 = -12345;
      header_DecQ.t1 = -12345;
      header_DecT.t1 = -12345;
      header_DecL.user4 = incnum;
      header_DecQ.user4 = incnum;
      header_DecT.user4 = incnum;
      header_DecL.user5 = data_DecL[maxN];
      header_DecQ.user5 = data_DecQ[maxN];
      header_DecT.user5 = data_DecT[maxN];
      string sl_name, sq_name, st_name, nsl_name, nsq_name, nst_name;
      sl_name =
          myFormatStringByFun((char *)"rf-%02.f-%03.f.sl.sac", incnum, win_len);
      sq_name =
          myFormatStringByFun((char *)"rf-%02.f-%03.f.sq.sac", incnum, win_len);
      st_name =
          myFormatStringByFun((char *)"rf-%02.f-%03.f.st.sac", incnum, win_len);
      nsl_name = myFormatStringByFun((char *)"rf-%02.f-%03.f.nsl.sac", incnum,
                                     win_len);
      nsq_name = myFormatStringByFun((char *)"rf-%02.f-%03.f.nsq.sac", incnum,
                                     win_len);
      nst_name = myFormatStringByFun((char *)"rf-%02.f-%03.f.nst.sac", incnum,
                                     win_len);

      write_sac(sl_name.data(), header_DecL, data_DecL);
      write_sac(sq_name.data(), header_DecQ, data_DecQ);
      write_sac(st_name.data(), header_DecT, data_DecT);
      write_sac(nsl_name.data(), header_DecL, data_DecLN);
      write_sac(nsq_name.data(), header_DecQ, data_DecQN);
      write_sac(nst_name.data(), header_DecT, data_DecTN);

    } // LOOP win_len STOP
  }   // LOOP incinum STOP

  end = clock(); //结束时间
  cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s"
       << endl; //输出时间（单位：ｓ）
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++functions in the main
// func+++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/*---------------------------------------------------------------------------*/
void fold(float z[], float r[], long lth, float *zz, float *rr, float *zr)
/* computes sums of squares and cross products
 * parameters of routine
 * float     z[];    input; z-array
 * float     r[];    input; r-array
 * long       lth;    input; length of traces
 * float       *zz;    output; sum of z squares
 * float       *rr;    output; sum of r squares
 * float       *zr;    output; sum of r*z products
 */
{
  /* local variables */
  long i; /* counter */
  /* executable code */
  *zz = *rr = *zr = 0.0;
  for (i = 0; i < lth; i++) {
    *zz += *z * *z;
    *rr += *r * *r;
    *zr += (*z++) * (*r++);
  } /*endfor*/
} /* end of fold */

void get_spikefilter(float *start, int trclth, float ff[], float delta,
                     float reg)
// float in_data[], float delta, float lowdw, float hiwdw,

/* computes spiking filter
 * par 1:    trace
 * par 2:    lo time bound
 * par 3:    hi time bound
 */
{
  /* local variables */
  // void    *trc;          /* input trace */
  // void    *new;          /* new created trace */
  // float     lowdw, hiwdw;  /* time window */
  // float    reg=1;           /* control parameter */
  // // long     loidx, hiidx;  /* sample window */
  // float     t0pos;         /* position of spike in sec */
  // float     lthfac;        /* length factor */
  // //int      trclth;        /* length of window in samples */
  // int      outlth;         length of output trace in samples
  // //float   *f, *ac, *ccr; /* pointers to data arrays */
  // float   *mov;          /* moving pointer */
  // //float   *start;        /* start of sample window */
  // float   tmp;           /* scratch */
  // int      d, i;          /* counter */
  // int      t0;            /* position of spike */
  // float reg=3;
  //-------determin mass center--------
  int t0 = 2;
  int outlth = trclth * 0.5;
  float tmp = 0.;
  float t0pos = 0.;
  float *mov;
  int i = 0, d;
  if (t0 == 2) { /* "center of mass" */
    tmp = 0.;
    t0pos = 0.;
    for (i = 0; i < trclth; i++) {
      tmp += (float)i * Abs(start[i]);
      t0pos += Abs(start[i]);
    } /*endfor*/
    t0 = Nint(tmp / t0pos);
  }                   /*endif*/
  else if (t0 == 1) { /* set t0 at maximum */
    tmp = *start;
    t0 = 0;
    for (i = 1; i < trclth; i++)
      if (fabs(start[i]) > tmp) {
        tmp = fabs(start[i]);
        t0 = i;
      } /*endif*/
  } else if (t0 == 3) {
    t0 = trclth - 1;
  }
  t0pos = t0 * delta;
  float f0 = t0;
  if ((f0 - int(f0)) >= delta) {
    t0 = int(f0) + 1;
  } else {
    t0 = int(f0);
  }

  t0pos = t0 * delta;
  // printf( "spiking %f:%d\n", t0pos,t0 );
  float *ac = new float[trclth];
  float *f = new float[outlth];
  for (int i = 0; i < trclth; i++) {
    ac[i] = 0;
  }
  for (int i = 0; i < outlth; i++) {
    f[i] = 0;
  }
  float *ccr;
  ccr = ac + outlth;
  /*  compute autocorrelation */
  mov = ac - 1;
  for (d = 0; d < outlth; d++) { /* shift counter */
    *(++mov) = 0.;
    for (i = d; i < trclth; i++)
      /* sample counter */
      *mov += start[i] * start[i - d];
  } /*endfor*/
  *ac *= reg + 1.0;
  /* compute correlation of start[0..trclth-1] with spike at pos t0 */
  if (t0 < 0)
    t0 = 0;
  if (t0 >= outlth)
    t0 = outlth - 1;
  mov = ccr;
  for (i = 0; i <= t0; i++) {
    *mov++ = start[t0 - i];
  }
  for (i = t0 + 1; i < outlth; i++) {
    *mov++ = 0.;
  }
  mt_levinson(ac, ccr, f, outlth);
  for (int ii = 0; ii < outlth; ii++) {
    ff[ii] = *f;
    f++;
  }
} /* end of get_spikefilter */

void mt_levinson(float *r, float *g, float *f, int m)
/* Levinson algorithm for symmetric toeplitz-matrix "r[0..m-1]".
 *
 * parameters of routine
 * float     r[];     input; symmetric toeplitz-matrix r[0..m-1]
 * float     g[];     input; correlation of input with desired trace
 * float     f[];     output; computed inverse filter
 * int        m;       input; length of all traces
 * int     *status; output; return status
 */
{
  /* local variables */
  float *a, *b; /* pointer to workspace */
  int i, j, ii; /* counter */
  float gn, z1, z2;
  /* executable code */
  // a = (float *)2L*(long)m;//
  a = new float[2 * (long)m];
  // if  (Severe(status))  return;
  b = a + m;
  /* create offset-1 arrays */
  r--;
  g--;
  f--;
  a--;
  b--;
  f[1] = g[1] / r[1];
  a[1] = r[2] / r[1];
  for (i = 2; i <= m; i++) {
    gn = r[1];
    z1 = (i == m) ? 0. : r[i + 1];
    z2 = g[i];
    for (j = 2; j <= i; j++) {
      gn -= r[j] * a[j - 1];
      z1 -= r[j] * a[i - j + 1];
      z2 -= r[j] * f[i - j + 1];
    } /*endfor 300*/
    a[i] = z1 / gn;
    f[i] = z2 / gn;
    ii = i - 1;
    for (j = 1; j <= ii; j++) {
      b[j] = a[j] - a[i] * a[ii - j + 1];
      f[j] -= f[i] * a[ii - j + 1];
    } /*endfor 400*/
    for (j = 1; j <= ii; j++) {
      a[j] = b[j];
    }
  } /*endfor 200*/
} /* end of mt_levinson */

void mt_fold(long la, float a[], long lb, float b[], float c[])
/* folds traces a[0..la-1] and b[0..lb-1] to trace c[0..la+lb-2]
 *
 * parameters of routine
 * long       la;       input; length of trace a
 * float       a[];      input; trace a
 * long       lb;       input; length of trace b
 * float       b[];      input; trace b
 * float       c[];      output; result trace (length la+lb-1)
 */
{
  /* local variables */
  long i, j, k, lc; /* counters */
  /* executable code */
  lc = la + lb - 1;
  /* clear output array */
  for (i = 0; i < lc; c[i++] = 0.0) {
  }
  /* compute convolution */
  for (i = 0; i < la; i++) {
    for (j = 0; j < lb; j++) {
      k = i + j;
      c[k] += a[i] * b[j];
    } /*endfor*/
      // sy_sharecpu();
  }   /*endfor*/
} /* end of mt_fold */

float sc_polar2(SAMPLE zarr_ori[], SAMPLE rarr_ori[], float delta,
                float mdir_t1, float mdir_t2)
/* computes main direction of polarisation in 2-dimensional case.
 * Result is the angle "angle" to the z-direction (in degrees)
 *
 * the coherence matrix
 *
 *          -               -        -        -
 *         |  <Z*Z>   <Z*R>  |      |  z    c  |
 *   C  =  |                 |  =:  |          |
 *         |  <Z*R>   <R*R>  |      |  c    r  |
 *          -               -        -        -
 *
 * is diagonalised.  The direction of the eigenvector with the
 * largest eigenvalue is the desired direction
 *
 * Let
 *   k := (z - r) / c
 * and
 *   alpha[1] := 1/2 * (k + sqrt(k*k+4))
 *   alpha[2] := 1/2 * (k - sqrt(k*k+4))
 *
 * then the eigenvalues lambda[i] and eigenvectors v[i] (i = 1,2)
 * are given by
 *
 *   lambda[i] = r + c*alpha[i]
 *   v[i] = ( alpha[i], 1 )
 *
 *
 * parameters of routine
 * SAMPLE     zarr[];  input; samples of z-trace
 * SAMPLE     rarr[];  input; samples of r-trace
 * long       lth;     input; length of traces in samples
 * REAL       *angle;  output; computed angle
 */
{
  /* local variables */
  REAL z, c, r, k; /* variables from above formula */
  REAL alpha[2];   /* see above */
  REAL lambda[2];  /* eigenvalues  of coherence matrix */
  REAL *zarr, *rarr;
  float angle = 0.;
  int mdir_t1_idx = int(mdir_t1 / delta);
  int mdir_t2_idx = int(mdir_t2 / delta);
  long lth = int((mdir_t2 - mdir_t1) / delta + 1);
  /* executable code */
  zarr = zarr_ori + mdir_t1_idx;
  rarr = rarr_ori + mdir_t1_idx;
  sc_trccorr(zarr, rarr, lth, &z, &r, &c);
  if (fabs(c) < SHC_EPSILON) {
    angle = 0.0;
  } else {
    k = (z - r) / c;
    alpha[0] = 0.5 * (k + sqrt(k * k + 4));
    alpha[1] = 0.5 * (k - sqrt(k * k + 4));
    lambda[0] = r + c * alpha[0];
    lambda[1] = r + c * alpha[1];
    if (fabs(lambda[0]) > fabs(lambda[1])) {
      angle = atan(1.0 / alpha[0]);
    } else {
      angle = atan(1.0 / alpha[1]);
    } /*endif*/
  }   /*endif*/

  angle *= SHC_RAD_TO_DEG;
  return angle;

} /* end of sc_polar2 */

void mt_rot2(SAMPLE *n, SAMPLE *e, long lth, REAL angle, SAMPLE r_out[],
             SAMPLE t_out[])
// void mt_rot2( SAMPLE *z, SAMPLE *r, long lth, REAL angle, SAMPLE *l, SAMPLE
// *q )
/* 2-dimensional rotation
 *
 * parameter of routine
 * SAMPLE    *n, *e;      input; north & east component
 * long      lth;         input; length of input & output traces
 * REAL      angle;       input; rotation angle
 * SAMPLE    *r, *t;      output; rotated traces
 */
{
  /* local variables */
  REAL ann, ane, aen, aee; /* rotation matrix */
  SAMPLE *r, *t;
  long i; /* counter */
  /* executable code */
  angle *= (SHC_PI / 180.);
  ann = cos(angle);
  ane = -sin(angle);
  aen = -ane;
  aee = ann;
  cout << "I am here!!!" << endl;
  for (i = 0; i <= lth; i++) {
    //*r++ = ann * *n + ane * *e;
    //*t++ = aen * *n + aee * *e;
    r_out[i] = ann * *n + ane * *e;
    t_out[i] = aen * *n + aee * *e;
    e++;
    n++;
  } /*endfor*/

} /* end of mt_rot2 */

void sc_trccorr(REAL z[], REAL r[], long lth, REAL *zz, REAL *rr, REAL *zr)
/* computes sums of squares and cross products
 *
 * parameters of routine
 * SAMPLE     z[];    input; z-array
 * SAMPLE     r[];    input; r-array
 * long       lth;    input; length of traces
 * REAL       *zz;    output; sum of z squares
 * REAL       *rr;    output; sum of r squares
 * REAL       *zr;    output; sum of r*z products
 */
{
  /* local variables */
  long i; /* counter */
  /* executable code */
  *zz = *rr = *zr = 0.0;
  for (i = 0; i < lth; i++) {
    *zz += *z * *z;
    *rr += *r * *r;
    *zr += (*z++) * (*r++);
  } /*endfor*/
} /* end of sc_trccorr */

void mt_rot3(SAMPLE *z, SAMPLE *n, SAMPLE *e, long lth, REAL azim, REAL inci,
             int type, SAMPLE l[], SAMPLE q[], SAMPLE t[])

/* 3-dimensional rotation
 *
 * parameter of routine
 * SAMPLE    *z, *n, *e;  input; vertical, north & east component
 * long      lth;         input; length of input & output traces
 * REAL      azim, inci;  input; rotation angles
 * int       type;        input; type of rotation
 * SAMPLE    *l, *q, *t;  output; rotated traces
 * type 1 SHC_ROT_ZNE_TO_LQT
 * type 2 SHC_ROT_ZNE_TO_UVW
 * type 3 SHC_ROT_UVW_TO_ZNE
 */
{
  /* local variables */
  REAL azz, azn, aze; /* rotation matrix */
  REAL anz, ann, ane;
  REAL aez, aen, aee;
  REAL cosi, cosa, sini, sina;
  long i; /* counter */

  /* executable code */

  azim *= (SHC_PI / 180.);
  cosa = cos(azim);
  sina = sin(azim);
  if (type == 1) {
    inci *= (SHC_PI / 180.);
    cosi = cos(inci);
    sini = sin(inci);
  } /*endif*/

  if (type == 1) {
    azz = cosi;
    azn = -sini * cosa;
    aze = -sini * sina;
    anz = sini;
    ann = cosi * cosa;
    ane = cosi * sina;
    aez = 0.0;
    aen = sina;
    aee = -cosa;
  } else if (type == 2) {
#ifdef XXX
    if (sina != 0.0)
      sina = 1.0 / sina;
    if (cosa != 0.0)
      cosa = 1.0 / cosa;
    azz = cosa / 3.0;
    anz = cosa / 3.0;
    aez = cosa / 3.0;
    azn = 0.;
    ann = sina / SQRT3;
    aen = -sina / SQRT3;
    aze = -sina * 2.0 / 3.0;
    ane = sina / 3.0;
    aee = sina / 3.0;
#endif
    azz = SQRT1_3;
    anz = SQRT1_3;
    aez = SQRT1_3;
    azn = 0.0;
    ann = SQRT1_2;
    aen = -SQRT1_2;
    aze = -SQRT2_3;
    ane = 0.5 * SQRT2_3;
    aee = 0.5 * SQRT2_3;
  } else if (type == 3) {
#ifdef XXX
    azz = cosa;
    anz = cosa;
    aez = cosa;
    azn = 0.;
    ann = sina * SQRT3_2;
    aen = -sina * SQRT3_2;
    aze = -sina;
    ane = sina * 0.5;
    aee = sina * 0.5;
#endif
    azz = SQRT1_3;
    azn = SQRT1_3;
    aze = SQRT1_3;
    anz = 0.0;
    ann = SQRT1_2;
    ane = -SQRT1_2;
    aez = -SQRT2_3;
    aen = 0.5 * SQRT2_3;
    aee = 0.5 * SQRT2_3;
  } else {
    azz = azn = aze = anz = ann = ane = aez = aen = aee = 0.0;
  } /*endif*/

  for (i = 0; i < lth; i++) {
    *l++ = azz * *z + azn * *n + aze * *e;
    *q++ = anz * *z + ann * *n + ane * *e;
    *t++ = aez * *z + aen * *n + aee * *e;
    z++;
    n++;
    e++;
  } /*endfor*/

} /* end of mt_rot3 */

char *myFormatStringByFun(char *format, ...) {
  va_list list;
  // 1. 先获取格式化后字符串的长度
  va_start(list, format);
  int size = vsnprintf(NULL, 0, format, list);
  va_end(list);
  if (size <= 0) {
    return NULL;
  }
  size++;

  // 2. 复位va_list，将格式化字符串写入到buf
  va_start(list, format);
  char *buf = (char *)malloc(size);
  vsnprintf(buf, size, format, list);
  va_end(list);
  return buf;
}