 #ifndef DYNAMICMEMORY_H
 #define DYNAMICMEMORY_H

#include <stdlib.h>
#include<math.h>

#ifdef __cplusplus
 extern "C" {
#endif

void* newTable1D(int len, int sizeOneElement);

void** newTable2D(int firstDim, int secondDim, int sizeOneElement);

void*** newTable3D(int firstDim, int secondDim, int thirdDim, int sizeOneElement);

void deleteTable1D(void* tablePtr);

void deleteTable2D(void** tablePtr);

void deleteTable3D(void*** tablePtr);

#ifdef __cplusplus
 }
#endif

 #endif
