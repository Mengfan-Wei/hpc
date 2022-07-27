/*
 *  This example illustrates how to create and close a group. 
 *  It is used in the HDF5 Tutorial.
 */

#include "hdf5.h"
#define FILE "group.h5"

int main() {

   hid_t       file_id, group_id, dataset_id, dataspace_id;  /* identifiers */
   hsize_t     dims[3];
   herr_t      status;

   /* Create a new file using default properties. */
   file_id = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

   /* Create a group named "/MyGroup" in the file. */
   group_id = H5Gcreate2(file_id, "/MyGroup", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   /* Create the data space for the dataset. */
   dims[0] = 5; 
   dims[1] = 6;
   dims[2] = 4;
   dataspace_id = H5Screate_simple(3, dims, NULL);

   /* Create a dataset named "/MyDataset" in the group. */
   dataset_id = H5Dcreate2(file_id, "/MyGroup/MyDataset", H5T_STD_I32BE, dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   /* End access to the dataset and release resources used by it. */
   status = H5Dclose(dataset_id);

   /* Close the group. */
   status = H5Gclose(group_id);

   /* Terminate access to the file. */
   status = H5Fclose(file_id);
}
