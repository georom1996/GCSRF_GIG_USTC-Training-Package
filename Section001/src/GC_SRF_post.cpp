/*------------------------------------------------------------------
----   Intially coded by Zhang Zhou @GIG on Jul., 7th, 2021   ----
----   Post Processing in GC_SRF strategy.  ----
-------------------------------------------------------------------*/

#include "sacio.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>

/*-------------------------------------------------- some CONST values
 * -----------------------------------*/
#define SAC_MAX 500 // Maximum of SAC format files.
#define LINE_LENGTH_MAX                                                        \
  200 // Maximum of characters of every line in "para.file".
char *myFormatStringByFun(char *format, ...);
float gama(float *A, float *B, int Length);
/*------------------------------------- get sign of float type number "x"
 * --------------------------------*/
float sign(float x) {
  if (x >= 0.)
    return 1.;
  else
    return -1.;
}
using namespace std;
/*---------------------------------------------- MAIN PROGRAM
 * --------------------------------------------*/
int main(int argc, char *argv[]) {
  clock_t start, end; // clock_t
  int i, j, k, shift_index, count = 1, npts, line = 0, tmp_index;
  float *data, *sum, delta, Nth_root, shift_time[SAC_MAX], weight[SAC_MAX],
      tmp1, tmp2, peak;
  char buff[1024], sac_out[200], s1[200], s2[200], *sac_name[SAC_MAX], norm[16];
  char *sacname, *min_name;
  float amp_0[16][20], energy_0[16][20];
  FILE *fin;
  SACHEAD hd;
  int incnum = 0;
  int incnum_good = 0;
  int ii, jj;
  ofstream outfile01, outfile02, outfile03, outfile04, outfile05, outfile06;
  outfile01.open("all.Amp_0.xyz");
  outfile02.open("all.Energy_0.xyz");
  outfile03.open("Min.lst");
  outfile04.open("all.CCC.xyz");
  outfile05.open("SelectedRed.CC.lst");
  outfile06.open("OptimalSRF.lst");

  //------------------------------LOOP BLOCK------------------------------------

  start = clock(); // Time recording
  int rank, comm_size;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  int N = 16;
  unsigned int box = N / comm_size;
  float win_len;
  float win_len_good = 0;

  for (win_len = 5 * (rank + 1); win_len <= 100;
       win_len = win_len + comm_size * 5) {
    int ii_minamp0 = 0, jj_minamp0 = 0;
    float min_amp = 100;
    int min_inci = 100;
    cout << win_len << "    " << rank << "    " << comm_size << endl;
    for (incnum = 0; incnum <= 60; incnum = incnum + 4) {
      float energy_tmp = 0;
      int idx_count = 0;

      ii = floor(incnum / 4);
      jj = floor(win_len / 5 - 1);
      sacname = myFormatStringByFun(
          (char *)"rf.%02d.%03.f.nsl.sac.100.cutrv.rm.rv", incnum, win_len);
      data = read_sac(sacname, &hd);
      delta = hd.delta;
      npts = hd.npts;
      int zeropoint = floor(100 / delta);
      amp_0[ii][jj] = data[zeropoint];

      if (fabs(amp_0[ii][jj]) <= fabs(min_amp)) {
        min_amp = amp_0[ii][jj];
        min_inci = incnum;
        min_name = sacname;
      }

      for (int t_index = floor((100 - 1) / delta);
           t_index <= floor((100 + 1) / delta); t_index++) {
        energy_tmp = energy_tmp + data[t_index] * data[t_index];
        idx_count++;
      }
      energy_0[ii][jj] = energy_tmp / idx_count;
      outfile01 << incnum << " " << win_len << " " << amp_0[ii][jj] << endl;
      outfile02 << incnum << " " << win_len << " " << energy_0[ii][jj] << endl;

    } // LOOP incinum STOP
    // cout<<"# mininci: "<<min_inci<<" min_amp: "<<min_amp<<" minlen:
    // "<<win_len<<"fname : "<<sacname<<endl;
    outfile03 << min_inci << " " << win_len << " " << min_amp << " " << min_name
              << endl;
    sac_name[jj] = min_name;
  } // LOOP win_len

  outfile01.close();
  outfile02.close();
  outfile03.close();

  // Get Ref_SRF=sum(Min_SRF)/N_win
  Nth_root = 1;
  line = jj;
  for (int j = 0; j <= line; j++) {
    shift_time[j] = 0;
    weight[j] = 1;
    peak = 0.;
    int norm_flag = 1;
    // cout<<"sac name:"<<sac_name[j]<<endl;
    data = read_sac(sac_name[j], &hd);
    delta = hd.delta;
    npts = hd.npts;
    if (j == 0) {
      sum = (float *)malloc(sizeof(npts) * npts);
      for (i = 0; i < npts; i++)
        sum[i] = 0.;
    }
    if (isnan(hd.depmax) != 1) {
      for (k = 0; k < npts; k++) {
        if (fabs(data[k]) > peak)
          peak = fabs(data[k]);
      }
      shift_index = (int)(shift_time[j] / delta);
      for (i = 0; i < npts; i++) {
        tmp_index = i + shift_index;
        if (tmp_index >= 0 && tmp_index < npts) {
          if (norm_flag == 1) {
            data[tmp_index] /= peak;
            tmp1 = sign(data[tmp_index]) *
                   pow(fabs(data[tmp_index]), 1. / Nth_root) * weight[j];
          } else {
            tmp1 = sign(data[tmp_index]) *
                   pow(fabs(data[tmp_index]), 1. / Nth_root) * weight[j];
          }
          sum[i] += tmp1 / line;
        } else
          sum[i] += 0.;
      }
    } else
      continue;
  }
  for (i = 0; i < npts; i++) {
    tmp2 = sign(sum[i]) * pow(fabs(sum[i]), Nth_root);
    sum[i] = tmp2;
    // cout<<sum[i]<<endl;
  }

  // CrossCorelation Coeff
  float maxcc = 0.0;
  char *optsrfname;
  for (win_len = 5; win_len <= 100; win_len = win_len + 5) {
    for (incnum = 0; incnum <= 60; incnum = incnum + 4) {
      float ccc = 0.0;
      ii = floor(incnum / 4);
      jj = floor(win_len / 5 - 1);
      sacname = myFormatStringByFun(
          (char *)"rf.%02d.%03.f.nsl.sac.100.cutrv.rm.rv", incnum, win_len);
      data = read_sac(sacname, &hd);
      delta = hd.delta;
      npts = hd.npts;
      ccc = gama(sum, data, npts);
      // cout<<"ccc :"<<ccc<<endl;
      outfile04 << incnum << " " << win_len << " " << ccc << endl;
      if (strcmp(sacname, sac_name[jj]) == 0) {
        if (ccc >= maxcc) {
          maxcc = ccc;
          optsrfname = sacname;
          incnum_good = incnum;
          win_len_good = win_len;
        }
        outfile05 << sacname << " " << incnum << " " << win_len << " " << ccc
                  << endl;
      }
    } // LOOP incinum STOP
  }   // LOOP win_len

  outfile06 << optsrfname << " " << incnum_good << " " << win_len_good << " "
            << maxcc << endl;
  cout << "GC_SRF done! Optimal SRF is : " << optsrfname << endl;
  outfile04.close();
  outfile05.close();
  outfile06.close();
  MPI_Finalize();

  end = clock(); // Time record end
  cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s"
       << endl; // print time on screen
}

//-------------------------------------------------------------------------------------------------
//
float gama(float *A, float *B, int Length) {

  float sumA = 0.0, sumB = 0.0, aveA = 0.0, aveB = 0.0;

  //
  for (int i = 0; i < Length; i++) {
    sumA = sumA + A[i];
    sumB = sumB + B[i];
  }

  //
  aveA = sumA / float(Length);
  aveB = sumB / float(Length);

  //
  double R1(0), R2(0), R3(0);
  for (int i = 0; i < Length; i++) {
    R1 += (A[i] - aveA) * (B[i] - aveB);
    R2 += pow((A[i] - aveA), 2);
    R3 += pow((B[i] - aveB), 2);
  }

  return (R1 / sqrt(R2 * R3));
}

//-------------------------------------------------------------------------------------------------

char *myFormatStringByFun(char *format, ...) {
  va_list list;
  // 1. 先获取格式化后字符串的长�?
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
