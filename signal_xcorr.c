/*
Copyright (c) 2019, Julien Chavanton
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL Julien Chavanton BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <math.h>
#include "stdint.h"
#include "kissfft/tools/kiss_fftr.h"
#include "kissfft/_kiss_fft_guts.h"

const int sampling_hz=8000;
int debug = 0;

FILE* open_file (char *fn) {
	printf("[%s][%s]\n", __FUNCTION__, fn);
	FILE *ref = fopen(fn, "r");
	if (!ref) return NULL;
	return ref;
}

/*
	https://dsp.stackexchange.com/questions/736/how-do-i-implement-cross-correlation-to-prove-two-audio-files-are-similar
	Cross-correlation and convolution are closely related. In short, to do convolution with FFTs, you

	zero-pad the input signals (add zeros to the end so that at least half of the wave is "blank")
	take the FFT of both signals
	multiply the results together (element-wise multiplication)
	do the inverse FFT
	conv(a, b) = ifft(fft(a_and_zeros) * fft(b_and_zeros))

	You need to do the zero-padding because the FFT method is actually circular cross-correlation,
	meaning the signal wraps around at the ends. So you add enough zeros to get rid of the overlap, to simulate a signal that is zero out to infinity.

	To get cross-correlation instead of convolution, you either need to time-reverse one of the signals before doing the FFT,
	or take the complex conjugate of one of the signals after the FFT:

	corr(a, b) = ifft(fft(a_and_zeros) * fft(b_and_zeros[reversed]))
	corr(a, b) = ifft(fft(a_and_zeros) * conj(fft(b_and_zeros))) << using complex conjugate
*/

float do_kissfft_xcorr (int nfft, char *fn, char *fn_ref, int *lag) {
	kiss_fftr_cfg fft_cfg, ffti_cfg;
	kiss_fft_scalar * rbuf_deg; // degraded
	kiss_fft_cpx * cbuf_deg;
	kiss_fft_scalar * rbuf_ref; // reference
	kiss_fft_cpx * cbuf_ref;
	kiss_fft_scalar * rbuf_res; // result
	kiss_fft_cpx * cbuf_res;
	int is_inverse = 0;

	fft_cfg = kiss_fftr_alloc(nfft , 0, NULL, NULL);
	ffti_cfg = kiss_fftr_alloc(nfft, 1, NULL, NULL);

	rbuf_deg = (kiss_fft_scalar*)malloc(sizeof(kiss_fft_scalar) * nfft * 2);
	cbuf_deg = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * (nfft/2+1) * 2);
	memset(rbuf_deg, sizeof(kiss_fft_scalar)*nfft*2, 0);

	rbuf_ref = (kiss_fft_scalar*)malloc(sizeof(kiss_fft_scalar) * nfft * 2);
	cbuf_ref = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * (nfft/2+1) * 2);
	memset(rbuf_ref, sizeof(kiss_fft_scalar)*nfft*2, 0);

	rbuf_res = (kiss_fft_scalar*)malloc(sizeof(kiss_fft_scalar) * nfft * 2);
	cbuf_res = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * (nfft/2+1) * 2);
	memset(rbuf_res, sizeof(kiss_fft_scalar)*nfft*2, 0);

	float fz_bin = sampling_hz / (float)(nfft);
	FILE *f = open_file(fn_ref);
	if (!f) return 0.0;
	FILE *f_ref = open_file(fn);
	if (!f) return 0.0;
	int y=0;
	int st = 1500;
		// read the degraded samples
		int x = 0;
		int16_t sample;
		int32_t sum = 0;
		while (fread(&sample, sizeof(int16_t), 1, f) > 0 ) {
			rbuf_deg[x] = (float)sample;
			sum += (sample*sample);
			x++;
			if (x == nfft) break;
		}
		if (x != nfft) {
			printf("read error[%d][%d]\n", x, nfft);
			while(x!=nfft) {
				rbuf_deg[x] = 0.0;
				x++;
			}
		}
		// read the reference samples
		x =0;
		y = 0;
		while (fread(&sample, sizeof(int16_t), 1, f_ref) > 0 ) {
			rbuf_ref[x] = (float)sample;
			x++;
			if (x == nfft) break;
		}
		if (x != nfft) {
			printf("read error[%d][%d]\n", x, nfft);
			while(x!=nfft) {
				rbuf_ref[x] = 0.0;
				x++;
			}
		}
		kiss_fftr(fft_cfg, rbuf_deg, cbuf_deg);
		kiss_fftr(fft_cfg, rbuf_ref, cbuf_ref);

		for (x=0;x<nfft;x++) {
			cbuf_ref[x].i = -cbuf_ref[x].i; // complex conjugate
			C_MUL(cbuf_res[x], cbuf_deg[x], cbuf_ref[x]);
			//printf ("kssfft [%d]xcorr[%f]\n", x, cbuf_res[x].r);
		}

		kiss_fftri(ffti_cfg, cbuf_res, rbuf_res);

		float max_corr = 0;
		for (x=0;x<nfft;x++) {
			if (max_corr < rbuf_res[x]) {
				max_corr = rbuf_res[x];
				*lag = x;
			}
		}
		printf("max xcorr[%.0f] sample[%d]\n", max_corr, *lag);

	free(ffti_cfg);
	free(fft_cfg);
	free(rbuf_ref);
	free(cbuf_ref);
	return max_corr;
}



int main(void) {
	printf("\nresults:\n\n");
	char *reference_fn = "files/reference.raw";
	char *degraded_fn = "files/degraded_lag100ms.raw";
	int lag = 0;
	int sampling_rate = 8000;
	float auto_corr, auto_corr_1, auto_corr_2, cross_corr;
//	int nfft = 8192;
//	int nfft = 16384;
	int nfft = 32768;
	auto_corr = do_kissfft_xcorr(nfft, reference_fn, reference_fn, &lag);
	cross_corr = do_kissfft_xcorr(nfft, degraded_fn, reference_fn, &lag);
	printf("coefficient correlation[%dms][%.0f/%.0f]: %.4f\n\n", (nfft-lag)/(sampling_rate/1000), cross_corr, auto_corr, cross_corr/auto_corr);

//	nfft = 32768;
	nfft = 16384;
	nfft = 8192;
	reference_fn = "files/kamailio_ref.raw";
	degraded_fn = "files/kamailio_1.raw";
	auto_corr_1 = do_kissfft_xcorr(nfft, degraded_fn, degraded_fn, &lag);
	auto_corr_2 = do_kissfft_xcorr(nfft, reference_fn, reference_fn, &lag);
	auto_corr = sqrt(auto_corr_1 * auto_corr_2);
	cross_corr = do_kissfft_xcorr(nfft, degraded_fn, reference_fn, &lag);
	printf("coefficient correlation[%dms][%.0f/%.0f]: %.4f\n\n", (nfft-lag)/(sampling_rate/1000), cross_corr, auto_corr, cross_corr/auto_corr);

	nfft = 8192;
	reference_fn = "files/10_1.raw";
	degraded_fn = "files/10_2.raw";
	auto_corr_1 = do_kissfft_xcorr(nfft, degraded_fn, degraded_fn, &lag);
	auto_corr_2 = do_kissfft_xcorr(nfft, reference_fn, reference_fn, &lag);
	auto_corr = sqrt(auto_corr_1 * auto_corr_2);
	cross_corr = do_kissfft_xcorr(nfft, degraded_fn, reference_fn, &lag);
	printf("coefficient correlation[%dms][%.0f/%.0f]: %.4f\n\n", (nfft-lag)/(sampling_rate/1000), cross_corr, auto_corr, cross_corr/auto_corr);

	nfft = 8192;
	reference_fn = "files/10_1.raw";
	degraded_fn = "files/dd_2.raw";
	auto_corr_1 = do_kissfft_xcorr(nfft, degraded_fn, degraded_fn, &lag);
	auto_corr_2 = do_kissfft_xcorr(nfft, reference_fn, reference_fn, &lag);
	auto_corr = sqrt(auto_corr_1 * auto_corr_2);
	cross_corr = do_kissfft_xcorr(nfft, degraded_fn, reference_fn, &lag);
	printf("coefficient correlation[%dms][%.0f/%.0f]: %.4f\n\n", (nfft-lag)/(sampling_rate/1000), cross_corr, auto_corr, cross_corr/auto_corr);

	nfft = 8192;
	nfft = 16384;
	nfft = 32768;
//	nfft = 65536;
	sampling_rate = 16000;
	reference_fn = "files/dd_1.raw";
	degraded_fn = "files/dd_2.raw";
	auto_corr_1 = do_kissfft_xcorr(nfft, degraded_fn, degraded_fn, &lag);
	auto_corr_2 = do_kissfft_xcorr(nfft, reference_fn, reference_fn, &lag);
	auto_corr = sqrt(auto_corr_1 * auto_corr_2);
	cross_corr = do_kissfft_xcorr(nfft, degraded_fn, reference_fn, &lag);
	printf("coefficient correlation[%dms][%.0f/%.0f]: %.4f\n\n", (nfft-lag)/(sampling_rate/1000), cross_corr, auto_corr, cross_corr/auto_corr);

	nfft = 16384;
	nfft = 32768;
	nfft = 65536;
	sampling_rate = 16000;
	reference_fn = "files/123_1.raw";
	degraded_fn = "files/123_2.raw";
	auto_corr_1 = do_kissfft_xcorr(nfft, degraded_fn, degraded_fn, &lag);
	auto_corr_2 = do_kissfft_xcorr(nfft, reference_fn, reference_fn, &lag);
	auto_corr = sqrt(auto_corr_1 * auto_corr_2);
	cross_corr = do_kissfft_xcorr(nfft, degraded_fn, reference_fn, &lag);
	printf("coefficient correlation[%dms][%.0f/%.0f]: %.4f\n\n", (nfft-lag)/(sampling_rate/1000), cross_corr, auto_corr, cross_corr/auto_corr);
	return 1;
}



