#include "HMM.h"
#include <assert.h>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "../Utilities.h"

int underflowProtection = 1; 
int haploid_underflowProtection = 0;
double underflowProtectionMultiplicator = 20;

