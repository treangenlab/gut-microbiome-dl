#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<ctype.h>
#include<math.h>

#include "devtest.h"

int k_to_ull_stride(int k) {
    return (int)ceil(k/32.);
}
