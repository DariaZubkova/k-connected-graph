#pragma once
#include <vector>
#include<iostream>

typedef struct node {
	int i;
	int j;
}node;

typedef struct edge {
	int id;
	int i;
	int j;
}edge;

class Matrix {
public:
	Matrix(int numLine, int numColumn);
	Matrix() {}
	void PrintMatrix();
	int currentNumEdge();
	int getElem(int i, int j);
	void setElem(int i, int j, int value);
	void addVal(int i, int j, int value);
	void difVal(int i, int j, int value);
	int getNumLine();
	int getNumColumn();
	void setNumLine(int num);
	void setNumColumn(int num);
	~Matrix();
	std::vector<int> arr;
private:
	int numLine;
	int numColumn;
};