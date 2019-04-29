#include <stdio.h>
#include <math.h>
// #define NRANSI
#include "stdint.h"
#include "kissfft/tools/kiss_fftr.h"
#include "kissfft/_kiss_fft_guts.h"
//#include "kiss_xcorr/include/kiss_xcorr.h"

const int sampling_hz=8000;
const int samples=4096; // must be muliple of 2 // 32768;
int N, NN, NN2;
int debug=0;

FILE* open_file (char *fn) {
	printf("[%s][%s]\n", __FUNCTION__, fn);
	FILE *ref = fopen(fn, "r");
	return ref;
}


// KISS FFT
static void read_raw_image(const char *fname, size_t N, kiss_fft_scalar *dst) {
	FILE *fp = fopen(fname, "rb");
	fread(dst, sizeof(kiss_fft_scalar), N, fp);
	fclose(fp);
}

int do_kissfft (char *fn) {
	kiss_fftr_cfg st;
	kiss_fft_scalar * rbuf; // float
	kiss_fft_cpx * cbuf;
	int is_inverse = 0;
	int nfft = 1024;

	rbuf = (kiss_fft_scalar*)malloc(sizeof(kiss_fft_scalar) * nfft );
	cbuf = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * (nfft/2+1) );
	st = kiss_fftr_alloc(nfft ,is_inverse ,0,0);

	float fz_bin = sampling_hz / (float)(nfft);
	FILE *f_reference = open_file(fn);
	int y=0;
	{
		int x =0;
		int16_t sample;
		while (fread(&sample, sizeof(int16_t), 1, f_reference) > 0 ) {
			rbuf[x] = (float)sample/32767.0;
			x++;
			if (x == nfft) break;
		}
		if (x != nfft) {
			printf("not enough samples in [%s]\n", fn);
			return 1;
		}
		kiss_fftr(st, rbuf, cbuf);

		int i;
		int bin = 0;
		float max_e = 0.0f, max_e_f = 0.0f;
		double norm_bin_mag, amplitude;
		for (i=0;i<(nfft/2);i++) {
			bin ++;
			if(cbuf[i].r > 0)
			if (((float)cbuf[i].r) > max_e) {
				max_e = (float)cbuf[i].r;
				max_e_f = bin*fz_bin;
				/* get normalized bin magnitude */
				// https://groups.google.com/forum/#!topic/comp.dsp/cZsS1ftN5oI
				norm_bin_mag = 2. * sqrt(cbuf[i].r*cbuf[i].r + cbuf[i].i*cbuf[i].i) / (float)nfft;
				norm_bin_mag = 2. * sqrt(cbuf[i].r*cbuf[i].r ) / (float)nfft;
				/* convert to dB value */
				amplitude = 20. * log10(norm_bin_mag);
			}
			norm_bin_mag = 2. * sqrt(cbuf[i].r*cbuf[i].r + cbuf[i].i*cbuf[i].i) / (float)nfft;
			amplitude = 20. * log10(norm_bin_mag);
			printf("d[%s] bin_fz[%f] mag[%f]dB[%f]\n", __FUNCTION__, bin*fz_bin, norm_bin_mag, amplitude);
		}
		printf("[%s] sampling_fz[%d] bin_fz[%f] max energy found in Fz[%f]\n", __FUNCTION__, sampling_hz, fz_bin, max_e_f);
	}

	free(st);
	free(rbuf);
	free(cbuf);
	return 1;
}

int test (){
	int isign;
	float *data;
	char *test_fn = "files/440.raw";
	do_kissfft(test_fn);
}

int main(void) {
	printf("\nresults:\n\n");
	test();
	return 1;
}
#undef NRANSI
