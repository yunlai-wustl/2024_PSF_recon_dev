#include <iostream>
#include <sys/stat.h>

#include "ImageArray.h"

template class ImageArray<int>;
template class ImageArray<float>;
//template class ImageArray<bool>;

template <typename T>
ImageArray<T>::ImageArray(){
	_dim_x = 0;
	_dim_y = 0;
	_dim_z = 0;
	_num_elem = (size_t)0;
	_image = NULL;
}

template <typename T>
ImageArray<T>::ImageArray(PET_geometry g){
	_image = NULL;
	Setup(g);
}

template <typename T>
ImageArray<T>::ImageArray(int dim1, int dim2, int dim3){
	_image = NULL;
	Setup(dim1, dim2, dim3);
}

template <typename T>
ImageArray<T>::ImageArray(ImageArray<T> & img, int x_downSampled_size, int y_downSampled_size, int z_downSampled_size) {
    // Setup the dimensions of the downsampled image
    _dim_x = img._dim_x / x_downSampled_size;
    _dim_y = img._dim_y / y_downSampled_size;
    _dim_z = img._dim_z / z_downSampled_size;
    _num_elem = _dim_x * _dim_y * _dim_z;
	_image = NULL;

    // Allocate memory for the downsampled image
    _AllocateData();

    // Iterate over the dimensions of the downsampled image
    for (int k = 0; k < _dim_z; ++k) {
        for (int j = 0; j < _dim_y; ++j) {
            for (int i = 0; i < _dim_x; ++i) {
                T sum = 0;
                int count = 0;

                // Iterate over the corresponding region in the original image
                for (int z = k * z_downSampled_size; z <min( img._dim_z, (k + 1) * z_downSampled_size); ++z) {
                    for (int y = j * y_downSampled_size; y < min(img._dim_y,(j + 1) * y_downSampled_size); ++y) {
                        for (int x = i * x_downSampled_size; x < min(img._dim_x, (i + 1) * x_downSampled_size); ++x) {
                            sum += img(x, y, z);
                            count++;
                        }
                    }
                }
                // Assign the average value to the downsampled image
                _image[(k * _dim_y + j) * _dim_x + i] = sum / count;
            }
        }
    }
}

template <typename T>
ImageArray<T>::~ImageArray(){
	_DeAllocateData();
}

template <typename T>
void
ImageArray<T>::Setup(PET_geometry g){
	_dim_x = g.NUM_X;
	_dim_y = g.NUM_Y;
	_dim_z = g.NUM_Z;
	_num_elem = _dim_x*_dim_y*_dim_z;

	_AllocateData();
}


template <typename T>
void
ImageArray<T>::Setup(int dim1, int dim2, int dim3){
	_dim_x = dim1;
	_dim_y = dim2;
	_dim_z = dim3;
	_num_elem = _dim_x*_dim_y*_dim_z;

	_AllocateData();
}

template <typename T>
void
ImageArray<T>::Reset(){
	_dim_x = 0;
	_dim_y = 0;
	_dim_z = 0;
	_num_elem = 0;

	_DeAllocateData();
}



template <typename T>
void
ImageArray<T>::_AllocateData(){
	if (_image != NULL){
		printf("Allocation failed, the image array is already allocated.\n");
		throw 1;
	}

	if ((_image = new T[_num_elem]) == NULL){
		printf("Allocation failed, can NOT allocate image array.\n");
		throw 1;
	}

}

template <typename T>
void
ImageArray<T>::_DeAllocateData(){
	if (_image != NULL){
		delete[] _image;
	}
}


template <typename T>
void
ImageArray<T>::ReadFromMem(T*& image){
	if (_image == NULL){
		printf("The image data field is not allocated yet. Allocation must be done before reading from file.\n");
		throw 1;
	}

	//This implementation is NOT safe! Please improve in the future.
	int i;
	for (i = 0; i <static_cast<int>( _num_elem); i++){
		_image[i] = image[i];
	}

}


template <typename T>
void
ImageArray<T>::AddFromMem(T*& image){

	if (_image == NULL){
		printf("The image data field is not allocated yet. Allocation must be done before reading from file.\n");
		throw 1;
	}

	//This implementation is NOT safe! Please improve in the future.
	int i;
	for (i = 0; i <static_cast<int>( _num_elem); i++){
		_image[i] += image[i];
	}
}




template <typename T>
void
ImageArray<T>::ReadFromFile(string filename){
	
	if (_image == NULL){
		printf("The image data field is not allocated yet. Allocation must be done before reading from file.\n");
		throw 1;
	}
	
	FILE *fid = fopen(filename.c_str(), "rbS");
	if (fid == NULL){
		printf("Could not open file %s . This file may not exist or is opened by another program.\n",filename.c_str());
		throw 1;
	}
	else
	{
		printf("......Reading %d elements from file %s\n", static_cast<int>(_num_elem), filename.c_str());
	}

	size_t count = fread(_image, sizeof(T), _num_elem, fid);
	fclose(fid);

	if (count != _num_elem){
		printf("Image file size does not match expected FOV size.\n");
		throw 1;
	}
}

template <typename T>
void
ImageArray<T>::WriteToFile(string filename){
	
	if (_image == NULL){
		printf("The image data field is not allocated yet.\n");
		throw 1;
	}
	
	FILE *fid = fopen(filename.c_str(), "wbS");
	if (fid == NULL){
		printf("Could not open file %s . This file may not exist or is opened by another program.\n");
		throw 1;
	}
	else
	{
		printf("......Writting %d elements from file %s\n", _num_elem, filename.c_str());
	}
	size_t count = fwrite(_image, sizeof(T), _num_elem, fid);
	fclose(fid);

	if (count != _num_elem){
		printf("The files written does not have the correct size.\n");
		throw 1;
	}
}


template <typename T>
void
ImageArray<T>::CopyFromImageArray(const ImageArray<T>& source){

	if (_image == NULL || _dim_x != source._dim_x || _dim_y != source._dim_y || _dim_z != source._dim_z){
		printf("ReAllocation needs to be done because this iamge array is not allocated or the source image array has different dimensions.\n");
		Reset();
		Setup(source._dim_x, source._dim_y, source._dim_z);
	}

	memcpy(_image, source._image, _num_elem*sizeof(T));

}


template <typename T>
void
ImageArray<T>::AddFromImageArray(const ImageArray<T>& source){

	if (_image == NULL || _dim_x != source._dim_x || _dim_y != source._dim_y || _dim_z != source._dim_z){
		printf("ReAllocation needs to be done because this iamge array is not allocated or the source image array has different dimensions.\n");
		Reset();
		Setup(source._dim_x, source._dim_y, source._dim_z);
	}

	for (int i = 0; i <static_cast<int>( _num_elem); i++)
	{
		_image[i] += source._image[i];
	}
}


template <typename T>
void
ImageArray<T>::MultiplyBy(const ImageArray<T>& source){

	if (_image == NULL || _dim_x != source._dim_x || _dim_y != source._dim_y || _dim_z != source._dim_z){
		printf("ReAllocation needs to be done because this iamge array is not allocated or the source image array has different dimensions.\n");
		Reset();
		Setup(source._dim_x, source._dim_y, source._dim_z);
	}

	for (int i = 0; i <static_cast<int>( _num_elem); i++)
	{
		_image[i] *= source._image[i];
	}
}

template <typename T>
void
ImageArray<T>::DivideBy(const ImageArray<T>& source){

	if (_image == NULL || _dim_x != source._dim_x || _dim_y != source._dim_y || _dim_z != source._dim_z){
		printf("ReAllocation needs to be done because this iamge array is not allocated or the source image array has different dimensions.\n");
		Reset();
		Setup(source._dim_x, source._dim_y, source._dim_z);
	}

	for (int i = 0; i <static_cast<int>( _num_elem); i++)
	{
		if (source._image[i] == 0){
			_image[i] = 0;
		}
		else{
			_image[i] /= source._image[i];
		}
	}
}

template <typename T>
void
ImageArray<T>::DivideBy(const ImageArray<T>& source, const ImageArray<T>& nonzero_mask){

	if (_image == NULL || _dim_x != source._dim_x || _dim_y != source._dim_y || _dim_z != source._dim_z){
		printf("ReAllocation needs to be done because this iamge array is not allocated or the source image array has different dimensions.\n");
		Reset();
		Setup(source._dim_x, source._dim_y, source._dim_z);
	}

	for (int i = 0; i <static_cast<int>( _num_elem); i++)
	{
		if (nonzero_mask[i] > 0.1){
			_image[i] /= source._image[i];
		}
		else{
			_image[i] = 0;
		}
	}
}

template <typename T>
void
ImageArray<T>::ScaledBy(const T& factor){

	int i;

	if (_image == NULL){
		printf("The image data field is not allocated yet.\n");
		throw 1;
	}
	else{
		for (i = 0; i <static_cast<int>( _num_elem); i++)
		{
			_image[i] *= factor;
		}
	}
}

template <typename T>
void
ImageArray<T>::SetValue(const T value){

	int i;

	if (_image == NULL){
		printf("Allocation needs to be done before setting the image value.\n");
		throw 1;
	}
	else{
		for (i = 0; i <static_cast<int>( _num_elem); i++)
		{
			_image[i] = value;
		}
	}

}

template <typename T>
T
ImageArray<T>::GetSum(){
	
	int i;
	T sum;

	sum = (T)0;

	if (_image == NULL){
		printf("Allocation needs to be done before setting the image value.\n");
		throw 1;
	}
	else{
		for (i = 0; i <static_cast<int>( _num_elem); i++){
			sum += _image[i];
		}
	}

	return sum;
}

template <typename T>
size_t
ImageArray<T>::GetSize(){

	if (_image == NULL){
		printf("Allocation needs to be done before setting the image value.\n");
		throw 1;
	}
	else{

	}

	return _num_elem;
}

template <typename T>
int 
ImageArray<T>::getDimX(){
	return _dim_x;
}

template <typename T>
int 
ImageArray<T>::getDimY(){
	return _dim_y;
}

template <typename T>
int 
ImageArray<T>::getDimZ(){
	return _dim_z;
}