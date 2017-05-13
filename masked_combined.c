#include "stdint.h"
#include "masked_combined.h"
#include "maths.h"
static uint8_t state[16][NUM_SHARES];
static uint8_t round_key[16][NUM_SHARES];

void Encrypt(uint8_t* output, uint8_t* input, uint8_t* key) {
    uint8_t rcon[10] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36};
    uint8_t round;
    uint8_t i, j, d;
    uint8_t tmp[4][NUM_SHARES];
        
    // Mask both the input and key
    MaskArray(state, input, 16);
    MaskArray(round_key, key, 16); 


    for(round = 0; round < 10; round++) {
        // Add Round Key Stage
        for(j = 0; j < NUM_SHARES; j++) {
            for(i = 0; i < 16; i++) {
                state[i][j] ^= round_key[i][j];
            }
        }
    
        CombinedSbox(state);
        ShiftRowsMixColumns(state, round);

        // Key Schedule
        round_key[0][0] ^= rcon[round];
        DualSbox(tmp[0], tmp[1], round_key[13], round_key[14]);
        DualSbox(tmp[2], tmp[3], round_key[15], round_key[12]);
        for(d = 0; d < NUM_SHARES; d++) {
            round_key[0][d] ^= tmp[0][d];
            round_key[1][d] ^= tmp[1][d];
            round_key[2][d] ^= tmp[2][d];
            round_key[3][d] ^= tmp[3][d];

            round_key[4][d] ^= round_key[0][d];
            round_key[5][d] ^= round_key[1][d];
            round_key[6][d] ^= round_key[2][d];
            round_key[7][d] ^= round_key[3][d];
            round_key[8][d] ^= round_key[4][d];
            round_key[9][d] ^= round_key[5][d];
            round_key[10][d] ^= round_key[6][d];
            round_key[11][d] ^= round_key[7][d];
            round_key[12][d] ^= round_key[8][d];
            round_key[13][d] ^= round_key[9][d];
            round_key[14][d] ^= round_key[10][d];
            round_key[15][d] ^= round_key[11][d];
        }
    }
    // Final Add Round Key
    for(j = 0; j < NUM_SHARES; j++) {
        for(i = 0; i < 16; i++) {
            state[i][j] ^= round_key[i][j];
        }
    }

    // Unmask the state revealing the encrypted output
    UnMaskArray(output, state, 16);
}

/*
    Takes an array and produces a set of shares for each element
*/
void MaskArray(uint8_t y[][NUM_SHARES], uint8_t x[], uint8_t length) {
    uint8_t i,j;
    for(i = 0; i < length; i++) {
        y[i][0] = x[i];
        for(j = 1; j < NUM_SHARES; j++) {
            y[i][j] = getRand();
            y[i][0] = y[i][0] ^ y[i][j];
        }
    }
}

/*
    Unmasked an array.
*/
void UnMaskArray(uint8_t y[], uint8_t x[][NUM_SHARES], uint8_t length) {
    uint8_t i,j;
    for(i = 0; i < length; i++) {
        y[i] = x[i][0];
        for(j = 1; j < NUM_SHARES; j++) {
            y[i] ^= x[i][j];
        }
    }
}


/*
 * Combined Sbox function
 * Computes the sbox for all 16 bytes of the state in 'parallel'
 */
void CombinedSbox(uint8_t s[16][NUM_SHARES]) {
    uint8_t i,j;
    uint8_t a[16][NUM_SHARES];
    uint8_t w[16][NUM_SHARES];
    
    SecEvalCombined(w, s, fifth);
    SecEvalCombined(a, w, fifth);
    SecEvalCombined(w, a, fifth);

    SecMultCombined(a, s, w, snd);

    for(i = 0; i < 16; i++) {
        for(j = 0; j < NUM_SHARES; j++) {
            // We have merged the affine lookup table with the last linear
            // squaring in the extended addition chain.
            s[i][j] = l_affine_snd[a[i][j]];
        }
        if((NUM_SHARES & 1) == 0) {
            s[i][0] ^= 0x63;
        }
    }
}



void DualSbox(uint8_t y1[], uint8_t y2[], uint8_t x1[], uint8_t x2[]) {

    uint8_t j;
    uint8_t w[2][NUM_SHARES];
    uint8_t a[2][NUM_SHARES];

    SecEvalTwoQuadraticRand(w[0], w[1], x1, x2, fifth);
    SecEvalTwoQuadraticRand(a[0], a[1], w[0], w[1], fifth);
    SecEvalTwoQuadraticRand(w[0], w[1], a[0], a[1], fifth);

    SecTwoMultSmall(a[0], a[1], x1, x2, w[0], w[1], snd);

    for(j = 0; j < NUM_SHARES; j++) {
        y1[j] = l_affine_snd[a[0][j]];
        y2[j] = l_affine_snd[a[1][j]];
    }
    if((NUM_SHARES & 1) == 0) {
        y1[0] ^= 0x63;
        y2[0] ^= 0x63;
    }
}

void ShiftRowsMixColumns(uint8_t s[][NUM_SHARES], uint8_t round) {
    uint8_t temp, share;
    uint8_t i;
    uint8_t Tmp,Tm,t;
    uint8_t r;
    for(share = 0; share < NUM_SHARES; share++) {
        // Rotate first row 1 columns to left
        temp         = s[1][share];
        s[1][share]  = s[5][share];
        s[5][share]  = s[9][share];
        s[9][share]  = s[13][share];
        s[13][share] = temp;

        // Rotate second row 2 columns to left
        temp         = s[2][share];
        s[2][share]  = s[10][share];
        s[10][share] = temp;

        temp         = s[6][share];
        s[6][share]  = s[14][share];
        s[14][share] = temp;

        // Rotate third row 3 columns to left
        temp         = s[3][share];
        s[3][share]  = s[15][share];
        s[15][share] = s[11][share];
        s[11][share] = s[7][share];
        s[7][share]  = temp;
  
        if(round < 9) {
            for(i = 0; i < 16; i += 4) {
                r = i;
                t   = s[r][share];
                Tmp = s[r][share] ^ s[r+1][share] ^ s[r+2][share] ^ s[r+3][share] ;
                
                Tm  = s[r][share] ^ s[r+1][share];
                Tm = xtime(Tm);
                s[r+0][share] ^= Tm ^ Tmp;

                Tm  = s[r+1][share] ^ s[r+2][share];
                Tm = xtime(Tm);
                s[r+1][share] ^= Tm ^ Tmp;

                Tm  = s[r+2][share] ^ s[r+3][share];
                Tm = xtime(Tm);
                s[r+2][share] ^= Tm ^ Tmp;

                Tm  = s[r+3][share] ^ t;
                Tm = xtime(Tm);
                s[r+3][share] ^= Tm ^ Tmp;
            }
        }
    }
}

/*
 * SecEval Function
 * Based on the CPRR SecEval function, this first runs the input 
 * through the Common Shares function and then uses the first 
 * technique by Zhang et al. for Random Reduction by 50%
 */
void SecEvalCombined(uint8_t w[16][NUM_SHARES], uint8_t z[16][NUM_SHARES], const uint8_t h[]) {
    uint8_t i,j,k,r,s0,t0,t1;

    for(i = 0; i < (NUM_SHARES/2); i++) {
        r = getRand();
        for(j = 0; j < 16; j++) {
            // Common Shares
            z[j][(NUM_SHARES/2) + i] = (z[j][(NUM_SHARES/2) + i] ^ r) ^ z[j][i];
            z[j][i] = r;
    
            // The first section of the O(n) space complexity CPRR SecEval
            w[j][(NUM_SHARES/2) + i] = h[z[j][(NUM_SHARES/2) + i]];        
            w[j][i] = h[r];
        }
    }

        
    for(i = 0; i < NUM_SHARES; i++) {
        for(j = (i + 1); j < NUM_SHARES; j++) {
            s0 = getRand();
            t0 = h[s0] ^ h[z[0][i] ^ s0];
            t1 = h[(z[0][i] ^ s0) ^ z[0][j]] ^ h[z[0][j] ^ s0];
            w[0][i] ^= t0;
            w[0][j] ^= t1;
 
            // Checks if the function can re-use the values already calculated
            if((i < (NUM_SHARES/2)) && (j < (NUM_SHARES/2))) {
                for(k = 1; k < 16; k++) {
                    w[k][i] ^= t0;
                    w[k][j] ^= t1;
                }
            } else {
                for(k = 1; k < 16; k++) {
                    s0 = getRand();
                    w[k][i] ^= h[s0] ^ h[z[k][i] ^ s0];
                    w[k][j] ^= h[(z[k][i] ^ s0) ^ z[k][j]] ^ h[z[k][j] ^ s0];
                }
            }
        }
    }
}

/*
 *  SecMult for 16 bytes
 */
void SecMultCombined(uint8_t c[][NUM_SHARES], uint8_t a[][NUM_SHARES], uint8_t b[][NUM_SHARES], const uint8_t f[]) {
    uint8_t tmp0, k, i, j;
    // First zero out the output array
    for(i = 0; i < NUM_SHARES; i++) {
        for(k = 0; k < 16; k++) {
            c[k][i] = 0;
        }
    }
    // Multiply the masked a and b for each of the 16 bytes
    for(i = 0; i < NUM_SHARES; i++) {
        for(k = 0; k < 16; k++) {
            c[k][i] ^= gfMul(f[a[k][i]],b[k][i]);
            for(j = (i+1); j < NUM_SHARES; j++) {
                tmp0 = getRand(); 
                c[k][i] ^=  tmp0;
                c[k][j] ^=  (tmp0 ^ gfMul(f[a[k][i]],b[k][j])) ^ gfMul(f[a[k][j]],b[k][i]);
            }
        }
    }
}

void SecTwoMultSmall(uint8_t c0[], uint8_t c1[], uint8_t aa0[], uint8_t aa1[], uint8_t bb0[], uint8_t bb1[], const uint8_t f[]) {
    uint8_t tmp0, tmp1, i, j;
    uint8_t ran = 0;
    for(i = 0; i < NUM_SHARES; i++) {
        c0[i] = gfMul(f[aa0[i]],bb0[i]);
        c1[i] = gfMul(f[aa1[i]],bb1[i]);
    }
    for(i = 0; i < NUM_SHARES; i++) {
        for(j = (i+1); j < NUM_SHARES; j++) {
            tmp0 = getRand(); ran++;
            c0[i] = c0[i] ^ tmp0;
            c0[j] = c0[j] ^ (tmp0 ^ gfMul(f[aa0[i]],bb0[j])) ^ gfMul(f[aa0[j]],bb0[i]);

            tmp1 = getRand(); ran++;
            c1[i] = c1[i] ^ tmp1;
            c1[j] = c1[j] ^ (tmp1 ^ gfMul(f[aa1[i]],bb1[j])) ^ gfMul(f[aa1[j]],bb1[i]);
        }
    }
}

void SecEvalTwoQuadraticRand(uint8_t y0[], uint8_t y1[], uint8_t aa0[], uint8_t aa1[], const uint8_t h[]) {
    uint8_t i, j, s0, s1, t0, t1, r;
    uint8_t a0[NUM_SHARES];
    uint8_t a1[NUM_SHARES];
    for(i = 0; i < (NUM_SHARES/2); i++) {
        r = getRand();
        a0[i] = r;
        a1[i] = r;

        a0[(NUM_SHARES/2) + i] = (aa0[(NUM_SHARES/2) + i] ^ r) ^ aa0[i];
        a1[(NUM_SHARES/2) + i] = (aa1[(NUM_SHARES/2) + i] ^ r) ^ aa1[i];

        y0[i] = h[r];
        y0[(NUM_SHARES/2) + i] = h[a0[(NUM_SHARES/2) + i]];
        y1[i] = h[r];
        y1[(NUM_SHARES/2) + i] = h[a1[(NUM_SHARES/2) + i]];
    }

    for(i = 0; i < NUM_SHARES; i++) {
        for(j = (i + 1); j < NUM_SHARES; j++) {
            s0 = getRand();
            t0 = h[s0] ^ h[a0[i] ^ s0];
            t1 = h[(a0[i] ^ s0) ^ a0[j]] ^ h[a0[j] ^ s0];
            y0[i] ^= t0;
            y0[j] ^= t1;
 
            if((i < (NUM_SHARES/2)) && (j < (NUM_SHARES/2))) {
                y1[i] ^= t0;
                y1[j] ^= t1;
            } else {
                s1 = getRand();
                y1[i] ^= h[s1] ^ h[a1[i] ^ s1];
                y1[j] ^= h[(a1[i] ^ s1) ^ a1[j]] ^ h[a1[j] ^ s1];
            }
        }
    }
}



