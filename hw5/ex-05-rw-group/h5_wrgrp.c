/*
 *  This example illustrates how to create and close a group. 
 *  It is used in the HDF5 Tutorial.
 */

#include "stdio.h"
#include "stdlib.h"
#include "hdf5.h"
#define FILE "group.h5"

int main() {

   hid_t       file_id, group_id, dataset_id;  /* identifiers */
   hsize_t     dims[3];
   herr_t      status;
   int         i, j, k,  dset_data[5][6][4];
   int * vec1 = (int*)malloc(5*6*4*sizeof(int));


   file_id = H5Fopen(FILE, H5F_ACC_RDWR, H5P_DEFAULT);
   dataset_id = H5Dopen(file_id, "/MyGroup/MyDataset", H5P_DEFAULT);

   status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
   printf("read from h5 dset_data[3][1][2]:%2d\n", dset_data[3][1][2]);

   status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec1);
   printf("vec1(3,1,2):%2d\n", vec1[6*4*3+4*1+2]);

   for (j = 0; j < 5; j++) {
    	for (i = 0; i < 6; i++){
      		for (k = 0; k < 4; k++)
			dset_data[j][i][k] = 200 + 100*j + 10*i + k;
     }
   }
   free(vec1);
   status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);

   status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
   printf("dset_data[3][1][2]:%2d\n", dset_data[3][1][2]);

   status = H5Dclose(dataset_id);

   /* Terminate access to the file. */
   status = H5Fclose(file_id);
}
