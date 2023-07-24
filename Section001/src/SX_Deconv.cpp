//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
typedef float SAMPLE; /* data type of seismogram samples */
typedef float REAL;	  /* general floating point type */
#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include "stdio.h"
#include <stdlib.h>
#include "./sac_lpz.h"
#define PI_HERE 3.14159265358979323846
#define RAD_TO_DEG (180.0 / PI_HERE)
#define EPSILON 1.0e-25
#define Abs(x) ((x) < 0 ? -(x) : (x))
#define Nint(x) ((x) > 0 ? (int)((x) + 0.5) : (int)((x)-0.5))
extern "C"
{
#include "./sac_lpz.h"
}
using namespace std;
// functions
void fold(float z[], float r[], long lth, float *zz, float *rr, float *zr);
void cut_data(float in_data[], float delta, float begin_time, float end_time, float out_data[]);
void get_spikefilter(float *start, int trclth, float ff[], float delta, float reg);
void mt_levinson(float *r, float *g, float *f, int m);
void spiking(float *start, int trclth, float f[]);
void mt_fold(long la, float a[], long lb, float b[], float c[]);
void mt_rot2(SAMPLE *n, SAMPLE *e, long lth, REAL angle, SAMPLE r[], SAMPLE t[]);
float sc_polar2(SAMPLE zarr_ori[], SAMPLE rarr_ori[], float delta, float mdir_t1, float mdir_t2);
void sc_trccorr(REAL z[], REAL r[], long lth, REAL *zz, REAL *rr, REAL *zr);
void mt_rot3(SAMPLE *z, SAMPLE *n, SAMPLE *e, long lth, REAL azim, REAL inci, SAMPLE l[], SAMPLE q[], SAMPLE t[]);
int main(int argc, char *argv[])
{
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
	for (int i = 0; i < argc; i++)
	{
		option_str = argv[i];
		if (option_str == "-Z")
		{
			name_Z = argv[i + 1];
			i++;
			continue;
		}
		if (option_str == "-N")
		{
			name_N = argv[i + 1];
			i++;
			continue;
		}
		if (option_str == "-E")
		{
			name_E = argv[i + 1];
			i++;
			continue;
		}
		if (option_str == "-t1")
		{
			spike_t1_str = argv[i + 1];
			spike_t1 = atof(spike_t1_str.c_str());
			i++;
			continue;
		}
		if (option_str == "-t2")
		{
			spike_t2_str = argv[i + 1];
			spike_t2 = atof(spike_t2_str.c_str());
			i++;
			continue;
		}
		if (option_str == "-reg")
		{
			reg_str = argv[i + 1];
			regnum = atof(reg_str.c_str());
			i++;
			continue;
		}
		if (option_str == "-inc")
		{
			inc_str = argv[i + 1];
			incnum = atof(inc_str.c_str());
			i++;
			continue;
		}
		if (option_str == "-D")
		{
			detial_flag = "detail_info";
			i++;
			continue;
		}
	}
	if (name_Z == "" && name_N == "" && detial_flag == "")
	{
		cout << "######################################################################" << endl;
		cout << "#############   Program for the Spiking deconvolution    #############" << endl;
		cout << "############# Version 1.02 Date 2020/03/29 By GeoRom/zhou ############" << endl;
		cout << "#############      Contact :zhagnzhou1996@gmail.com     ##############" << endl;
		cout << "######################################################################" << endl;
		cout << "S deconv ERROR : Pleas check the options:" << endl;
		cout << "Right usage :" << endl;
		cout << "       SX_Deconv -Z [file name of Z] -N [file name of N] -E [file name of E]  [-t1 t_spike_begin] [-t2 t_spike_end] -reg [control fil] -inc [inci] -D [with internal files]" << endl;
		cout << "       e.g. : SX_Deconv -Z Z.sac -N N.sac -E E.sac -R R.sac -T T.sac -D -t1 t1 -t2 t2" << endl;
		cout << "       -t1 for the begin of the spike" << endl;
		cout << "       -t2 for the end of spike" << endl;
		cout << "       -D for output the filter file" << endl;
		cout << "Output files: rf.sl.sac rf.sq.sac rf.st.sac" << endl;
		return 1;
	}
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	float *data_Z;
	float delta_Z;
	long int npts_Z;
	SACHEAD header_Z;
	data_Z = read_sac(name_Z.c_str(), &header_Z);
	if (data_Z == 0)
	{
		cout << "!error when read Z sac data!" << endl;
		return -1;
	}
	delta_Z = header_Z.delta;
	npts_Z = header_Z.npts;
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	float *data_N;
	float delta_N;
	long int npts_N;
	SACHEAD header_N;
	data_N = read_sac(name_N.c_str(), &header_N);
	if (data_N == 0)
	{
		cout << "!error when read N sac data!" << endl;
		return -1;
	}
	delta_N = header_N.delta;
	npts_N = header_N.npts;
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	float *data_E;
	float delta_E;
	long int npts_E;
	SACHEAD header_E;
	data_E = read_sac(name_E.c_str(), &header_E);
	if (data_E == 0)
	{
		cout << "!error when read E sac data!" << endl;
		return -1;
	}
	delta_E = header_E.delta;
	npts_E = header_E.npts;
	SACHEAD header_L = header_Z, header_Q = header_Z, header_T = header_Z;
	float delta_L = header_L.delta, delta_Q = header_Q.delta, delta_T = header_T.delta;
	long npts_L = header_L.npts, npts_Q = header_Q.npts, npts_T = header_T.npts;
	strcpy(header_L.kcmpnm, "L-COM");
	strcpy(header_Q.kcmpnm, "Q-COM");
	strcpy(header_T.kcmpnm, "T-COM");
	float azimnum = header_Z.baz;
	float data_L[int(npts_L)], data_Q[int(npts_Q)], data_T[int(npts_T)];
	mt_rot3(data_Z, data_N, data_E, npts_Z, azimnum, incnum, data_L, data_Q, data_T);
	cout << "ZNE2LQT PROCESSING: BAZ=" << azimnum << " INCI=" << incnum << "    OK!" << endl;
	if (detial_flag != "")
	{
		write_sac("l-com.sac", header_L, data_L);
		write_sac("q-com.sac", header_Q, data_Q);
		write_sac("t-com.sac", header_T, data_T);
	}
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	SACHEAD header_fil = header_Z, header_DecL = header_L, header_DecQ = header_Q, header_DecT = header_T;
	float spike_b = 20;
	float spike_e = 80;
	if (spike_t1 != 0.0)
	{
		spike_b = spike_t1;
	}
	if (spike_t2 != 0.0)
	{
		spike_e = spike_t2;
	}
	int indexlo = 0;
	int indexhi = 0;
	indexlo = int(spike_b / header_Q.delta) + 1;
	indexhi = int(spike_e / header_Q.delta) + 1;

	int mv_points = int(indexhi - indexlo + 1);
	int npts_fil = int(mv_points / 2);
	float filter_spk[npts_fil];
	float *mov;
	mov = data_Q;
	get_spikefilter(mov + indexlo, mv_points, filter_spk, header_Q.delta, regnum);
	header_fil.npts = npts_fil;
	header_fil.b = 0;
	strcpy(header_fil.kcmpnm, "SPK_FIL");
	if (detial_flag != "")
	{
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
	for (int i = indexlo; i < indexhi; i++)
	{
		if (max < data_DecQ[i])
		{
			max = data_DecQ[i];
			maxN = i;
		}
	}
	strcpy(header_DecL.kcmpnm, "DecL-OUT");
	strcpy(header_DecQ.kcmpnm, "DecQ-OUT");
	strcpy(header_DecT.kcmpnm, "DecT-OUT");
	header_DecL.npts = npts_L + header_fil.npts - 1;
	header_DecQ.npts = npts_Q + header_fil.npts - 1;
	header_DecT.npts = npts_T + header_fil.npts - 1;
	for (int i = 0; i < header_DecL.npts; i++)
	{
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
	write_sac("rf.sl.sac", header_DecL, data_DecL);
	write_sac("rf.sq.sac", header_DecQ, data_DecQ);
	write_sac("rf.st.sac", header_DecT, data_DecT);
	write_sac("rf.nsl.sac", header_DecL, data_DecLN);
	write_sac("rf.nsq.sac", header_DecQ, data_DecQN);
	write_sac("rf.nst.sac", header_DecT, data_DecTN);
}
//++++++++++++++++++++++++++++++         Functions        +++++++++++++++++++++++++++++++++++++
void fold(float z[], float r[], long lth, float *zz, float *rr, float *zr)
{
	long i;
	*zz = *rr = *zr = 0.0;
	for (i = 0; i < lth; i++)
	{
		*zz += *z * *z;
		*rr += *r * *r;
		*zr += (*z++) * (*r++);
	}
}

void get_spikefilter(float *start, int trclth, float ff[], float delta, float reg)
{
	int t0 = 2;
	int outlth = trclth * 0.5;
	float tmp = 0.;
	float t0pos = 0.;
	float *mov;
	int i = 0, d;
	if (t0 == 2)
	{
		tmp = 0.;
		t0pos = 0.;
		for (i = 0; i < trclth; i++)
		{
			tmp += (float)i * Abs(start[i]);
			t0pos += Abs(start[i]);
		}
		t0 = Nint(tmp / t0pos);
	}
	else if (t0 == 1)
	{
		tmp = *start;
		t0 = 0;
		for (i = 1; i < trclth; i++)
			if (fabs(start[i]) > tmp)
			{
				tmp = fabs(start[i]);
				t0 = i;
			} /*endif*/
	}
	else if (t0 == 3)
	{
		t0 = trclth - 1;
	}
	t0pos = t0 * delta;
	float f0 = t0;
	if ((f0 - int(f0)) >= delta)
	{
		t0 = int(f0) + 1;
	}
	else
	{
		t0 = int(f0);
	}
	t0pos = t0 * delta;
	float *ac = new float[trclth];
	float *f = new float[outlth];
	for (int i = 0; i < trclth; i++)
	{
		ac[i] = 0;
	}
	for (int i = 0; i < outlth; i++)
	{
		f[i] = 0;
	}
	float *ccr;
	ccr = ac + outlth;

	mov = ac - 1;
	for (d = 0; d < outlth; d++)
	{
		*(++mov) = 0.;
		for (i = d; i < trclth; i++)
			*mov += start[i] * start[i - d];
	}
	*ac *= reg + 1.0;

	if (t0 < 0)
		t0 = 0;
	if (t0 >= outlth)
		t0 = outlth - 1;
	mov = ccr;
	for (i = 0; i <= t0; i++)
	{
		*mov++ = start[t0 - i];
	}
	for (i = t0 + 1; i < outlth; i++)
	{
		*mov++ = 0.;
	}
	mt_levinson(ac, ccr, f, outlth);
	for (int ii = 0; ii < outlth; ii++)
	{
		ff[ii] = *f;
		f++;
	}
}

void mt_levinson(float *r, float *g, float *f, int m)
{
	float *a, *b;
	int i, j, ii;
	float gn, z1, z2;
	a = new float[2 * (long)m];
	b = a + m;
	r--;
	g--;
	f--;
	a--;
	b--;
	f[1] = g[1] / r[1];
	a[1] = r[2] / r[1];
	for (i = 2; i <= m; i++)
	{
		gn = r[1];
		z1 = (i == m) ? 0. : r[i + 1];
		z2 = g[i];
		for (j = 2; j <= i; j++)
		{
			gn -= r[j] * a[j - 1];
			z1 -= r[j] * a[i - j + 1];
			z2 -= r[j] * f[i - j + 1];
		}
		a[i] = z1 / gn;
		f[i] = z2 / gn;
		ii = i - 1;
		for (j = 1; j <= ii; j++)
		{
			b[j] = a[j] - a[i] * a[ii - j + 1];
			f[j] -= f[i] * a[ii - j + 1];
		}
		for (j = 1; j <= ii; j++)
		{
			a[j] = b[j];
		}
	}
}

void mt_fold(long la, float a[], long lb, float b[], float c[])
{
	long i, j, k, lc;
	lc = la + lb - 1;
	for (i = 0; i < la; i++)
	{
		for (j = 0; j < lb; j++)
		{
			k = i + j;
			c[k] += a[i] * b[j];
		}
	}
}

float sc_polar2(SAMPLE zarr_ori[], SAMPLE rarr_ori[], float delta, float mdir_t1, float mdir_t2)
{
	REAL z, c, r, k;
	REAL alpha[2];
	REAL lambda[2];
	REAL *zarr, *rarr;
	float angle = 0.;
	int mdir_t1_idx = int(mdir_t1 / delta);
	int mdir_t2_idx = int(mdir_t2 / delta);
	long lth = int((mdir_t2 - mdir_t1) / delta + 1);
	zarr = zarr_ori + mdir_t1_idx;
	rarr = rarr_ori + mdir_t1_idx;
	sc_trccorr(zarr, rarr, lth, &z, &r, &c);
	if (fabs(c) < EPSILON)
	{
		angle = 0.0;
	}
	else
	{
		k = (z - r) / c;
		alpha[0] = 0.5 * (k + sqrt(k * k + 4));
		alpha[1] = 0.5 * (k - sqrt(k * k + 4));
		lambda[0] = r + c * alpha[0];
		lambda[1] = r + c * alpha[1];
		if (fabs(lambda[0]) > fabs(lambda[1]))
		{
			angle = atan(1.0 / alpha[0]);
		}
		else
		{
			angle = atan(1.0 / alpha[1]);
		}
	}

	angle *= RAD_TO_DEG;
	return angle;
}

void mt_rot2(SAMPLE *n, SAMPLE *e, long lth, REAL angle, SAMPLE r_out[], SAMPLE t_out[])
{
	REAL ann, ane, aen, aee;
	SAMPLE *r, *t;
	long i;
	angle *= (PI_HERE / 180.);
	ann = cos(angle);
	ane = -sin(angle);
	aen = -ane;
	aee = ann;
	for (i = 0; i <= lth; i++)
	{
		r_out[i] = ann * *n + ane * *e;
		t_out[i] = aen * *n + aee * *e;
		e++;
		n++;
	}
}

void sc_trccorr(REAL z[], REAL r[], long lth, REAL *zz, REAL *rr, REAL *zr)
{
	long i;
	*zz = *rr = *zr = 0.0;
	for (i = 0; i < lth; i++)
	{
		*zz += *z * *z;
		*rr += *r * *r;
		*zr += (*z++) * (*r++);
	}
}

void mt_rot3(SAMPLE *z, SAMPLE *n, SAMPLE *e, long lth, REAL azim, REAL inci, SAMPLE l[], SAMPLE q[], SAMPLE t[])
{
	REAL azz, azn, aze;
	REAL anz, ann, ane;
	REAL aez, aen, aee;
	REAL cosi, cosa, sini, sina;
	long i;
	azim *= (PI_HERE / 180.);
	cosa = cos(azim);
	sina = sin(azim);
	inci *= (PI_HERE / 180.);
	cosi = cos(inci);
	sini = sin(inci);
	azz = cosi;
	azn = -sini * cosa;
	aze = -sini * sina;
	anz = sini;
	ann = cosi * cosa;
	ane = cosi * sina;
	aez = 0.0;
	aen = sina;
	aee = -cosa;
	for (i = 0; i < lth; i++)
	{
		*l++ = azz * *z + azn * *n + aze * *e;
		*q++ = anz * *z + ann * *n + ane * *e;
		*t++ = aez * *z + aen * *n + aee * *e;
		z++;
		n++;
		e++;
	}
}
