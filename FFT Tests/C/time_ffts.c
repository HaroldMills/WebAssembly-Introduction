/*

C program that times real FFTs of various sizes.

The functions real_fft and _reverse_bin_permute are translations to C of
the functions realFft and reverseBinPermute of the JavaScript file
https://github.com/HaroldMills/Vesper/blob/master/vesper/django/app/static/vesper/dft.js.
Those functions are in turn derived from functions RFFT.prototype.forward
and reverseBinPermute of the JavaScript file dsp.js
(https://github.com/corbanbrook/dsp.js) by Corban Brook, which are derived
from the split_radix_real_complex_fft and revbin_permute functions of the
C++ library FXT (www.jjj.de/fxt) by Jörg Arndt.

Harold Mills
January 2018

*/


#include <math.h>
#include <stdio.h>
#include <time.h>


void _time_ffts();
int _fft_cost(int e);
void _get_fft_test_input(double *output, int n, int f);
void real_fft(double *input, int n, double *output);
void _reverse_bin_permute(double *dest, double *source, int n);
void _test_real_fft();
void _test_reverse_bin_permute();


int main(int argc, char ** argv) {
    _time_ffts();
    // _test_reverse_bin_permute();
    // _test_real_fft();
}


void _time_ffts() {

    int freq = 1;
    int min_exp = 4;
    int max_exp = 16;
    double max_fft_cost = _fft_cost(max_exp);
    int num_trials = 10;
    int min_num_ffts_per_trial = 1000;

    for (int e = min_exp; e <= max_exp; ++e) {

        int fft_size = 1 << e;
        double input[fft_size];
        double output[fft_size];

        _get_fft_test_input(input, fft_size, freq);

        // Get trial size.
        double fft_cost_ratio = max_fft_cost / _fft_cost(e);
        int num_ffts_per_trial = min_num_ffts_per_trial * fft_cost_ratio;

        double min_us_per_fft = 1e20;

        for (int i = 0; i != num_trials; ++i) {

            clock_t begin = clock();

            for (int j = 0; j != num_ffts_per_trial; ++j)
                real_fft(input, fft_size, output);

            clock_t end = clock();
            double elapsed = (end - begin) / ((double) CLOCKS_PER_SEC);
            double us_per_fft = 1e6 * elapsed / num_ffts_per_trial;

            if (us_per_fft < min_us_per_fft)
                min_us_per_fft = us_per_fft;

            // printf(
            //   "    %d %d %ld %ld %lf %lf\n", i, num_ffts_per_trial, begin, end,
            //   elapsed, us_per_fft);

        }

        // printf("%d %d %lf\n", fft_size, num_ffts_per_trial, min_us_per_fft);
        printf("%lf\n", min_us_per_fft);

    }

}


int _fft_cost(int n) {
    return n * (1 << n);
}


void _get_fft_test_input(double *output, int n, int f) {
    double factor = 2 * M_PI * f / n;
    for (int i = 0; i != n; ++i)
        output[i] = cos(factor * i);
}


void real_fft(double *input, int n, double *output) {

    double *x = output;
    double TWO_PI = 2 * M_PI;
    int n2, n4, n8, nn;
    double t1, t2, t3, t4;
    int ix, id, i0, i1, i2, i3, i4, i5, i6, i7, i8;
    double st1, cc1, ss1, cc3, ss3;
    double e,  a, rval, ival, mag;

    _reverse_bin_permute(output, input, n);

    for (int ix = 0, id = 4; ix < n; id *= 4) {
        for (int i0 = ix; i0 < n; i0 += id) {
            //sumdiff(x[i0], x[i0+1]); // {a, b}  <--| {a+b, a-b}
            st1 = x[i0] - x[i0+1];
            x[i0] += x[i0+1];
            x[i0+1] = st1;
        }
        ix = 2*(id-1);
    }

    n2 = 2;
    nn = n >> 1;

    while((nn = nn >> 1)) {
        ix = 0;
        n2 = n2 << 1;
        id = n2 << 1;
        n4 = n2 >> 2;
        n8 = n2 >> 3;
        do {
            if(n4 != 1) {
                for(i0 = ix; i0 < n; i0 += id) {
                    i1 = i0;
                    i2 = i1 + n4;
                    i3 = i2 + n4;
                    i4 = i3 + n4;

                    //diffsum3_r(x[i3], x[i4], t1); // {a, b, s} <--| {a, b-a, a+b}
                    t1 = x[i3] + x[i4];
                    x[i4] -= x[i3];
                    //sumdiff3(x[i1], t1, x[i3]);   // {a, b, d} <--| {a+b, b, a-b}
                    x[i3] = x[i1] - t1;
                    x[i1] += t1;

                    i1 += n8;
                    i2 += n8;
                    i3 += n8;
                    i4 += n8;

                    //sumdiff(x[i3], x[i4], t1, t2); // {s, d}  <--| {a+b, a-b}
                    t1 = x[i3] + x[i4];
                    t2 = x[i3] - x[i4];

                    t1 = -t1 * M_SQRT1_2;
                    t2 *= M_SQRT1_2;

                    // sumdiff(t1, x[i2], x[i4], x[i3]); // {s, d}  <--| {a+b, a-b}
                    st1 = x[i2];
                    x[i4] = t1 + st1;
                    x[i3] = t1 - st1;

                    //sumdiff3(x[i1], t2, x[i2]); // {a, b, d} <--| {a+b, b, a-b}
                    x[i2] = x[i1] - t2;
                    x[i1] += t2;
                }
            } else {
                for(i0 = ix; i0 < n; i0 += id) {
                    i1 = i0;
                    i2 = i1 + n4;
                    i3 = i2 + n4;
                    i4 = i3 + n4;

                    //diffsum3_r(x[i3], x[i4], t1); // {a, b, s} <--| {a, b-a, a+b}
                    t1 = x[i3] + x[i4];
                    x[i4] -= x[i3];

                    //sumdiff3(x[i1], t1, x[i3]);   // {a, b, d} <--| {a+b, b, a-b}
                    x[i3] = x[i1] - t1;
                    x[i1] += t1;
                }
            }

            ix = (id << 1) - n2;
            id = id << 2;
        } while (ix < n);

        e = TWO_PI / n2;

        for (int j = 1; j < n8; j++) {
          a = j * e;
          ss1 = sin(a);
          cc1 = cos(a);

          //ss3 = sin(3*a); cc3 = cos(3*a);
          cc3 = 4*cc1*(cc1*cc1-0.75);
          ss3 = 4*ss1*(0.75-ss1*ss1);

          ix = 0; id = n2 << 1;
          do {
            for (i0 = ix; i0 < n; i0 += id) {
              i1 = i0 + j;
              i2 = i1 + n4;
              i3 = i2 + n4;
              i4 = i3 + n4;

              i5 = i0 + n4 - j;
              i6 = i5 + n4;
              i7 = i6 + n4;
              i8 = i7 + n4;

              //cmult(c, s, x, y, &u, &v)
              //cmult(cc1, ss1, x[i7], x[i3], t2, t1); // {u,v} <--| {x*c-y*s, x*s+y*c}
              t2 = x[i7]*cc1 - x[i3]*ss1;
              t1 = x[i7]*ss1 + x[i3]*cc1;

              //cmult(cc3, ss3, x[i8], x[i4], t4, t3);
              t4 = x[i8]*cc3 - x[i4]*ss3;
              t3 = x[i8]*ss3 + x[i4]*cc3;

              //sumdiff(t2, t4);   // {a, b} <--| {a+b, a-b}
              st1 = t2 - t4;
              t2 += t4;
              t4 = st1;

              //sumdiff(t2, x[i6], x[i8], x[i3]); // {s, d}  <--| {a+b, a-b}
              //st1 = x[i6]; x[i8] = t2 + st1; x[i3] = t2 - st1;
              x[i8] = t2 + x[i6];
              x[i3] = t2 - x[i6];

              //sumdiff_r(t1, t3); // {a, b} <--| {a+b, b-a}
              st1 = t3 - t1;
              t1 += t3;
              t3 = st1;

              //sumdiff(t3, x[i2], x[i4], x[i7]); // {s, d}  <--| {a+b, a-b}
              //st1 = x[i2]; x[i4] = t3 + st1; x[i7] = t3 - st1;
              x[i4] = t3 + x[i2];
              x[i7] = t3 - x[i2];

              //sumdiff3(x[i1], t1, x[i6]);   // {a, b, d} <--| {a+b, b, a-b}
              x[i6] = x[i1] - t1;
              x[i1] += t1;

              //diffsum3_r(t4, x[i5], x[i2]); // {a, b, s} <--| {a, b-a, a+b}
              x[i2] = t4 + x[i5];
              x[i5] -= t4;
            }

            ix = (id << 1) - n2;
            id = id << 2;

          } while (ix < n);
        }
    }

    // Scale output to have same norm as input.
    double f = 1 / sqrt(n);
    for (int i = 0; i < n; i++)
  	    x[i] *= f;

};


void _reverse_bin_permute(double *dest, double *source, int n) {

    int n2 = n >> 1;
    int nm1 = n - 1;
    int i = 1, r = 0, h;

    dest[0] = source[0];

    do {

        r += n2;

        dest[i] = source[r];
        dest[r] = source[i];

        i++;

        h = n2 << 1;
        while (h >>= 1, !((r ^= h) & h));

        if (r >= i) {

          dest[i]     = source[r];
          dest[r]     = source[i];

          dest[nm1 - i] = source[nm1 - r];
          dest[nm1 - r] = source[nm1 - i];

        }

        i++;

    } while (i < n2);

    dest[nm1] = source[nm1];

}


void _test_reverse_bin_permute() {

    int n = 8;
    double x[n];
    double y[n];

    for (int i = 0; i != n; ++i)
        x[i] = i;

    _reverse_bin_permute(y, x, n);

    for (int i = 0; i != n; ++i)
        printf("%d %lf %lf\n", i, x[i], y[i]);

}


void _test_real_fft() {

    int fft_size = 8;
    int freq = 1;
    double input[fft_size];
    double output[fft_size];

    _get_fft_test_input(input, fft_size, freq);

    real_fft(input, fft_size, output);

    for (int i = 0; i != fft_size; ++i)
        printf("%d %lf\n", i, output[i]);

}
