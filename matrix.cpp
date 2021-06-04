#include "matrix.h"

void Matrix::printMatrix() {
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
	if (numLine > 0)
		return numLine;
	return 0;
}

int Matrix::getNumColumn() {
	if (numColumn > 0)
		return numColumn;
	return 0;
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

void Matrix::addVal(int i, int j, int value) {
	arr[numColumn * i + j] += value;
}

void Matrix::difVal(int i, int j, int value) {
	arr[numColumn * i + j] -= value;
}

int Matrix::getElem(int i, int j) {
	if (!arr.empty()) {
		if (numColumn > 0) {
			if (arr[numColumn * i + j] > 0) {
				return arr[numColumn * i + j];
			}
		}
	}
	return 0;
}

Matrix::Matrix(int numLine = 0, int numColumn = 0) {
	this->numLine = numLine;
	this->numColumn = numColumn;
	std::vector<int> arr(numLine * numColumn, 0);
	this->arr = arr;
}

Matrix::~Matrix() {
	if (!arr.empty()) {
		this->arr.clear();
	}
}