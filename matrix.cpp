#include "matrix.h"

void Matrix::PrintMatrix() {
	for (int i = 0; i < numLine; i++) {
		for (int j = 0; j < numColumn; j++) {
			std::cout << arr[numColumn * i + j] << " ";
		}
		std::cout << std::endl;
	}
}

int Matrix::currentNumEdge() {
	int numEdge = 0;
	for (int i = 0; i < numLine - 1; i++) {
		for (int j = i + 1; j < numColumn; j++) {
			if (arr[numColumn * i + j] == 1)
				numEdge++;
		}
	}
	return numEdge;
}

int Matrix::getNumLine() {
	return numLine;
}

int Matrix::getNumColumn() {
	return numColumn;
}

void Matrix::setNumLine(int num) {
	numLine = num;
}

void Matrix::setNumColumn(int num) {
	numColumn = num;
}

void Matrix::setElem(int i, int j, int value) {
	arr[numColumn * i + j] = value;
}

int Matrix::getElem(int i, int j) {
	return arr[numColumn * i + j];
}

Matrix::Matrix(int numLine = 0, int numColumn = 0) {
	this->numLine = numLine;
	this->numColumn = numColumn;
	std::vector<int> arr(numLine * numColumn);
	this->arr = arr;
}

Matrix::~Matrix() {
	this->arr.clear();
}