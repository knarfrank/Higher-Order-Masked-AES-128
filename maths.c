#include "maths.h"
#include "stdint.h"
uint8_t gfMul(uint8_t a, uint8_t b) {
    int s = 0;
#if GFMUL_ASM    
    __asm__ __volatile__ (
        "ldrb r0, [%[table], %[a]]\n\t" // r0 = table[a]
        "ldrb r1, [%[table], %[b]]\n\t" // r1 = table[b]
        "add r0, r0, r1\n\t"            // r0 = table[a] + table[b]


        "add %[table], %[table], #255\n\t"
        "add %[table], %[table], #1\n\t"
        "ldrb r0, [%[table], r0]\n\t" // r0 = atable[r0]; 

        // Check if a or b == 0
        "mov r1, %[a]\n\t"
        "mul r1, r1, %[b]\n\t"
        "beq set_zero\n\t"
        
        "mov %[s], r0\n\t"
        "b gfmul_return\n\t"
        
        "set_zero:"
        "mov %[s], #0\n\t"    
        "b gfmul_return\n\t"

        "gfmul_return:\n\t"
        :[s]"+l" (s)
        :[a]"l" (a), [b]"l" (b), [table]"l" (table)
        :"r0", "r1", "memory"    
    );
#else
    s = table[a] + table[b]; 
	int q;
    /* Get the antilog */
	s = table[s+256];
	uint8_t z = 0;
    q = s;
    if(a == 0) {
		s = z;
	} else {
		s = q;
    }
	if(b == 0) {
		s = z; 
	} else {
		q = z;
	}
#endif
	return s;
}

uint8_t getRand(void) {
    // This should be replace with actual random number generator that isn't terrible
    return rand() % 256;
}
