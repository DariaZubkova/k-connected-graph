#include"graph.h"
#include"matrix.h"
#include <map>
#include <queue>
#include <utility>
#include <string>
#include <ctime>
#include <fstream>
#include <thread>
#include <mutex>
#include <chrono>
#include <functional>

std::mutex mtx;
std::vector<Graph> array_graph;

bool compareMatrix(std::vector<Graph>& arr_graph, Matrix& matrix) {
	int nL = matrix.getNumLine(), nC = matrix.getNumColumn();
	bool res = true;
	//matrix.PrintMatrix();
	//std::cout << "COMPARE:" << std::endl;
	if (arr_graph.empty()) {
		return false;
	}
	for (auto v : arr_graph) {
		res = true;
		//v.PrintMatrix();
		//std::cout << std::endl;
		Matrix m = v.get_Matrix();
		for (int i = 0; i < nL; i++) {
			for (int j = 0; j < nC; j++) {
				if (m.getElem(i, j) != matrix.getElem(i, j)) {
					res = false;
				}
			}
		}
		if (res) {
			return res;
		}
	}
	//std::cout << std::endl;
	return res;
}

void fun_fun2(Graph& graph, int* num, int k, std::vector<Graph>* arr_graph) { //ïåðåáîð ãðàôîâ
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	std::vector<int> dVertex;
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex.push_back(s);
	}
	for (int i = 0; i < numLine - 1; i++) { //numLine - 1
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k && dVertex[j] > k) {
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				dVertex[i]--;
				dVertex[j]--;
				graph.set_Matrix(matrix);
				//std::cout << "Delete " << i << " " << j << std::endl;
				//(*num)++;
				//PrintMatrix(arr, numLine, numColumn);
				//check
				res = graph.algorithmEven();
				if (res && !compareMatrix((*arr_graph), matrix)) {
					//std::cout << *num + 1 << std::endl;
					(*arr_graph).push_back(graph);
				}
				(*num) += 1;
				fun_fun2(graph, num, k, arr_graph);
				dVertex[i]++;
				dVertex[j]++;
				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
			}
		}
	}
	dVertex.clear();
}


void fun(Graph& graph, int* num, int k, std::vector<Graph>* arr_graph) { //ïåðåáîð ãðàôîâ
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	std::vector<int> dVertex;
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex.push_back(s);
	}
	for (int i = 4; i < 5; i++) { //numLine - 1
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k&& dVertex[j] > k) {
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				dVertex[i]--;
				dVertex[j]--;
				graph.set_Matrix(matrix);
				//std::cout << "Delete " << i << " " << j << std::endl;
				//(*num)++;
				//PrintMatrix(arr, numLine, numColumn);
				//check
				res = graph.algorithmEven();
				if (res && !compareMatrix((*arr_graph), matrix)) {
					//std::cout << *num + 1 << std::endl;
					(*arr_graph).push_back(graph);
				}
				(*num) += 1;
				fun_fun2(graph, num, k, arr_graph);
				dVertex[i]++;
				dVertex[j]++;
				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
			}
		}
	}
	dVertex.clear();
}

/*
class Thread_Graph
{
public:
	Thread_Graph() {}
	~Thread_Graph() {}

	void thread_graph(Graph& graph, int k) {
		Matrix matrix = graph.get_Matrix();
		int numLine = matrix.getNumLine();
		int numColumn = matrix.getNumColumn();
		bool res = false;
		int maxFlow;
		std::vector<std::thread> threads_graph;
		std::vector<int> dVertex;
		for (int i = 0; i < numLine - 1; i++) {
			for (int j = i + 1; j < numColumn; j++) {
				Graph gr = graph;
				Matrix m = gr.get_Matrix();
				for (int i = 0; i < numLine; i++) {
					int s = 0;
					for (int j = 0; j < numColumn; j++) {
						s += m.getElem(i, j);
					}
					dVertex.push_back(s);
				}
				if (m.getElem(i, j) == 1 && dVertex[i] > k&& dVertex[j] > k) {
					m.setElem(i, j, 0);
					if (m.getElem(j, i) == 1) {
						m.setElem(j, i, 0);
					}
					dVertex[i]--;
					dVertex[j]--;
					gr.set_Matrix(m);
					//std::cout << "Delete " << i << " " << j << std::endl;
					//(*num)++;
					//PrintMatrix(arr, numLine, numColumn);
					//check
					res = gr.algorithmEven();
					if (res && !compareMatrix(array, m)) {
						//std::cout << *num + 1 << std::endl;

						array.push_back(gr);

					}

					std::thread t(&Thread_Graph::fun_thread, this, gr, k);
					threads_graph.push_back(std::move(t));
					//threads_graph.push_back(std::thread(fun_thread, gr, k, arr_graph));
					//fun_thread(gr, k, arr_graph);

					dVertex[i]++;
					dVertex[j]++;
					m.setElem(i, j, 1);
					m.setElem(j, i, 1);
					gr.set_Matrix(matrix);
				}
				dVertex.clear();
			}
		}
		//dVertex.clear();
		//for (int i = 0; i < threads_graph.size(); i++) {
		//	threads_graph[i].detach();
		//}
	}

	std::vector<Graph> get_array() {
		return array;
	}

private:

	std::vector<Graph> list;
	std::recursive_mutex mutex_graph_class;
	std::vector<Graph> array;

	void fun_thread(Graph& graph, int k) { //ïåðåáîð ãðàôîâ
		Matrix matrix = graph.get_Matrix();
		int numLine = matrix.getNumLine();
		int numColumn = matrix.getNumColumn();
		bool res = false;
		int maxFlow;
		std::vector<int> dVertex;
		for (int i = 0; i < numLine; i++) {
			int s = 0;
			for (int j = 0; j < numColumn; j++) {
				s += matrix.getElem(i, j);
			}
			dVertex.push_back(s);
		}
		for (int i = 0; i < numLine - 1; i++) {
			for (int j = i + 1; j < numColumn; j++) {
				if (matrix.getElem(i, j) == 1 && dVertex[i] > k&& dVertex[j] > k) {
					matrix.setElem(i, j, 0);
					if (matrix.getElem(j, i) == 1) {
						matrix.setElem(j, i, 0);
					}
					dVertex[i]--;
					dVertex[j]--;
					graph.set_Matrix(matrix);
					//std::cout << "Delete " << i << " " << j << std::endl;
					//(*num)++;
					//PrintMatrix(arr, numLine, numColumn);
					//check
					res = graph.algorithmEven();
					if (res && !compareMatrix(array, matrix)) {
						//std::cout << *num + 1 << std::endl;
						mutex_graph_class.lock();
						array.push_back(graph);
						mutex_graph_class.unlock();
					}

					fun_thread(graph, k);

					dVertex[i]++;
					dVertex[j]++;
					matrix.setElem(i, j, 1);
					matrix.setElem(j, i, 1);
					graph.set_Matrix(matrix);
				}
			}
		}
		dVertex.clear();
	}

};*/

//std::function <void(Graph&, int, std::vector<Graph>*)> fun_thread = [](Graph& graph, int k, std::vector<Graph>* arr_graph)


void fun_thread2(Graph graph, int k) { //ïåðåáîð ãðàôîâ
	Graph gr = graph;
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	std::vector<int> dVertex(numLine);
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex[i] = s;
	}
	for (int i = 0; i < numLine - 1; i++) {
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k&& dVertex[j] > k) {
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				dVertex[i]--;
				dVertex[j]--;
				graph.set_Matrix(matrix);
				//std::cout << "Delete " << i << " " << j << std::endl;
				//(*num)++;
				//PrintMatrix(arr, numLine, numColumn);
				//check
				res = graph.algorithmEven();
				if (res && !compareMatrix(array_graph, matrix)) {
					//std::cout << *num + 1 << std::endl;
					//mutex_graph.lock();
					array_graph.push_back(graph);
					//mutex_graph.unlock();
				}

				fun_thread2(graph, k);

				dVertex[i]++;
				dVertex[j]++;
				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
			}
		}
	}
	dVertex.clear();
}


void fun_thread(Graph graph, int k) { //ïåðåáîð ãðàôîâ
	//std::cout << nn << std::endl;
	//std::cout << "ID THREAD = " << std::this_thread::get_id() << std::endl;
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	std::vector<int> dVertex(numLine);
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex[i] = s;
	}
	for (int i = 0; i < numLine - 1; i++) {
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k&& dVertex[j] > k) {
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				graph.set_Matrix(matrix);
				//(*num)++;
				res = graph.algorithmEven();
				if (res && !compareMatrix(array_graph, matrix)) {
					//std::cout << *num + 1 << std::endl;
					std::lock_guard<std::mutex> lock(mtx);
					//mutex_graph.lock();
					array_graph.push_back(graph);
					//mutex_graph.unlock();
				}

				fun_thread(graph, k);

				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
			}
		}
	}
	dVertex.clear();
}

void f2(Graph graph, int k) {
	std::cout << "ID THREAD = " << std::this_thread::get_id() << std::endl;
	Graph gr = graph;
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	std::vector<int> dVertex(numLine);
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex.push_back(s);
	}
	for (int i = 0; i < numLine - 1; i++) {
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k&& dVertex[j] > k) {
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				dVertex[i]--;
				dVertex[j]--;
				graph.set_Matrix(matrix);
				//std::cout << "Delete " << i << " " << j << std::endl;
				//(*num)++;
				//PrintMatrix(arr, numLine, numColumn);
				//check
				res = graph.algorithmEven();
				if (res && !compareMatrix(array_graph, matrix)) {
					//std::cout << *num + 1 << std::endl;
					//mutex_graph.lock();
					array_graph.push_back(graph);
					//mutex_graph.unlock();
				}

				fun_thread2(graph, k);

				dVertex[i]++;
				dVertex[j]++;
				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
			}
		}
	}
	dVertex.clear();
}

int algo(Graph graph) {
	int res = graph.algorithmEven();
	return res;
}

void f3(Graph graph, int k) {
	//std::cout << "ID THREAD = " << std::this_thread::get_id() << std::endl;
	if (graph.matrix.arr.empty()) {
		std::cout << "GOVNO =  " << std::this_thread::get_id() << std::endl;
		return;
	}
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	std::vector<int> dVertex(numLine);
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex[i] = s;
	}
	for (int i = 0; i < numLine - 1; i++) {
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k && dVertex[j] > k) {
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				//std::cout << "i = " << i << std::endl;
				//std::cout << "j = " << j << std::endl;
				graph.set_Matrix(matrix);
				//(*num)++;
				//res = graph.algorithmEven();
				//res = algo(graph);
				if (res && !compareMatrix(array_graph, matrix)) {
					//std::cout << *num + 1 << std::endl;
					std::lock_guard<std::mutex> lock(mtx);
					//mutex_graph.lock();
					array_graph.push_back(graph);
					//mutex_graph.unlock();
				}

				//std::cout << "Next = " << std::this_thread::get_id() << std::endl;
				f3(graph, k);

				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
				//std::cout << "End = " << i << " " << j << std::endl;
			}
		}
	}
	dVertex.clear();
}

void f4(Graph graph, int k) {
	//std::cout << "ID THREAD = " << std::this_thread::get_id() << std::endl;
	if (graph.matrix.arr.empty()) {
		std::cout << "GOVNO =  " << std::this_thread::get_id() << std::endl;
		return;
	}
	std::vector<std::thread> threads;
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	std::vector<int> dVertex(numLine);
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex[i] = s;
	}
	for (int i = 0; i < numLine - 1; i++) {
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k && dVertex[j] > k) {
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				//std::cout << "i = " << i << std::endl;
				//std::cout << "j = " << j << std::endl;
				graph.set_Matrix(matrix);
				//(*num)++;
				//res = graph.algorithmEven();
				//res = algo(graph);
				if (res && !compareMatrix(array_graph, matrix)) {
					//std::cout << *num + 1 << std::endl;
					std::lock_guard<std::mutex> lock(mtx);
					//mutex_graph.lock();
					array_graph.push_back(graph);
					//mutex_graph.unlock();
				}

				//std::cout << "Next = " << std::this_thread::get_id() << std::endl;
				//f(graph, k);
				std::thread thr(f3, graph, k);
				//threads.push_back(std::move(thr));
				threads.emplace_back(std::move(thr));

				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
				//std::cout << "End = " << i << " " << j << std::endl;
			}
		}
	}
	//dVertex.clear();
	for (auto& thr : threads) {
		thr.join();
	}
}

void f5(Graph graph, int k) {
	//std::cout << "ID THREAD = " << std::this_thread::get_id() << std::endl;
	if (graph.matrix.arr.empty()) {
		std::cout << "GOVNO =  " << std::this_thread::get_id() << std::endl;
		return;
	}
	std::vector<std::thread> threads;
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	std::vector<int> dVertex(numLine);
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex[i] = s;
	}
	for (int i = 0; i < numLine - 1; i++) {
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k && dVertex[j] > k) {
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				//std::cout << "i = " << i << std::endl;
				//std::cout << "j = " << j << std::endl;
				graph.set_Matrix(matrix);
				//(*num)++;
				//res = graph.algorithmEven();
				//res = algo(graph);
				if (res && !compareMatrix(array_graph, matrix)) {
					//std::cout << *num + 1 << std::endl;
					std::lock_guard<std::mutex> lock(mtx);
					//mutex_graph.lock();
					array_graph.push_back(graph);
					//mutex_graph.unlock();
				}

				//std::cout << "Next = " << std::this_thread::get_id() << std::endl;
				//f(graph, k);
				std::thread thr(f4, graph, k);
				//threads.push_back(std::move(thr));
				threads.emplace_back(std::move(thr));

				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
				//std::cout << "End = " << i << " " << j << std::endl;
			}
		}
	}
	//dVertex.clear();
	for (auto& thr : threads) {
		thr.join();
	}
}

void f1(Graph graph, int k) {
	//std::cout << "ID THREAD = " << std::this_thread::get_id() << std::endl;
	if (graph.matrix.arr.empty()) {
		std::cout << "GOVNO =  " << std::this_thread::get_id() << std::endl;
		return;
	}
	std::vector<std::thread> threads;
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	std::vector<int> dVertex(numLine);
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex[i] = s;
	}
	for (int i = 0; i < numLine - 1; i++) {
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k && dVertex[j] > k) {
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				//std::cout << "i = " << i << std::endl;
				//std::cout << "j = " << j << std::endl;
				graph.set_Matrix(matrix);
				//(*num)++;
				//res = graph.algorithmEven();
				//res = algo(graph);
				if (res && !compareMatrix(array_graph, matrix)) {
					//std::cout << *num + 1 << std::endl;
					std::lock_guard<std::mutex> lock(mtx);
					//mutex_graph.lock();
					array_graph.push_back(graph);
					//mutex_graph.unlock();
				}

				//std::cout << "Next = " << std::this_thread::get_id() << std::endl;
				//f(graph, k);
				std::thread thr(f5, graph, k);
				//threads.push_back(std::move(thr));
				threads.emplace_back(std::move(thr));

				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
				//std::cout << "End = " << i << " " << j << std::endl;
			}
		}
	}
	//dVertex.clear();
	for (auto& thr : threads) {
		thr.join();
	}
}


void f(Graph graph, int k) {
	std::cout << "ID THREAD = " << std::this_thread::get_id() << std::endl;
	//if (graph.matrix.arr.empty()) {
	//	std::cout << "GOVNO =  " << std::this_thread::get_id() << std::endl;
	//	return;
	//}
	std::vector<std::thread> threads;
	std::cout << "matrix get" << std::endl;
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	std::cout << "dVertex init" << std::endl;
	std::vector<int> dVertex(numLine);
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex[i] = s;
	}
	std::cout << "For i" << std::endl;
	for (int i = 0; i < numLine - 1; i++) {
		std::cout << "For j" << std::endl;
		for (int j = i + 1; j < numColumn; j++) {
			std::cout << "If" << std::endl;
			matrix.PrintMatrix();
			//for (int re = 0; re < numLine; re++) {
			//	std::cout << "dVertex[" << re << "] = " << dVertex[re] << std::endl;
			//}
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k && dVertex[j] > k) { // problem!
				std::cout << "sets" << std::endl;
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				//std::cout << "i = " << i << std::endl;
				//std::cout << "j = " << j << std::endl;
				std::cout << "set matrix" << std::endl;
				graph.set_Matrix(matrix);
				std::cout << "set matrix end" << std::endl;
				//(*num)++;
				//std::cout << "lock" << std::endl;
				//std::lock_guard<std::mutex> lock(mtx);
				std::cout << "Even" << std::endl;
				res = graph.algorithmEven();
				//res = algo(graph);
				if (res && !compareMatrix(array_graph, matrix)) {
					//std::cout << *num + 1 << std::endl;
					//std::lock_guard<std::mutex> lock(mtx);
					//mutex_graph.lock();
					std::cout << "push res" << std::endl;
					array_graph.push_back(graph);
					//mutex_graph.unlock();
				}

				//std::cout << "Next = " << std::this_thread::get_id() << std::endl;
				std::cout << "recursive" << std::endl;
				f(graph, k);

				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
				//std::cout << "End = " << i << " " << j << std::endl;
			}
		}
	}
	dVertex.clear();
}


void thread_graph(Graph& graph, int k) {
	std::cout << "ID THREAD MAIN" << std::endl;
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	std::vector<std::thread> threads;
	std::vector<int> dVertex(numLine);
	for (int i = 0; i < numLine - 1; i++) { //numLine - 1
		for (int j = i + 1; j < numColumn; j++) { //numColumn
			Graph gr(graph.get_numNode(), graph.get_k());
			Matrix m(numLine, numColumn);
			for (int u = 0; u < numLine; u++) {
				for (int l = 0; l < numColumn; l++) {
					m.setElem(u,l, matrix.getElem(u, l));
				}
			}
			//Matrix m = graph.get_Matrix();
			gr.set_Matrix(m);
			for (int h = 0; h < numLine; h++) {
				int s = 0;
				for (int l = 0; l < numColumn; l++) {
					s += matrix.getElem(h, l);
				}
				dVertex[h] = s;
			}
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k&& dVertex[j] > k) {
				matrix.setElem(i, j, 0);
				m.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
					m.setElem(j, i, 0);
				}
				graph.set_Matrix(matrix);
				gr.set_Matrix(m);
				//std::cout << "Delete " << i << " " << j << std::endl;
				//(*num)++;
				//PrintMatrix(arr, numLine, numColumn);
				//check
				//res = graph.algorithmEven();
				res = gr.algorithmEven();
				//res = algo(graph);
				if (res && !compareMatrix(array_graph, m)) { //matrix
					//std::cout << *num + 1 << std::endl;
					std::lock_guard<std::mutex> lock(mtx);
					array_graph.push_back(gr); //graph

				}

				std::thread thr(f, gr, k); //graph
				//threads.push_back(std::move(thr));
				threads.emplace_back(std::move(thr));
				//threads.emplace_back(fun_thread, gr, k);
				//std::thread t(&fun_thread, std::ref(gr), k);
				//t.detach();
				//threads_graph.push_back(std::move(t));
				//threads_graph.push_back(std::thread(&fun_thread, std::ref(gr), k, arr_graph));
				//fun_thread(gr, k, arr_graph);

				matrix.setElem(i, j, 1);
				m.setElem(i, j, 0);
				matrix.setElem(j, i, 1);
				m.setElem(j, i, 0);
				graph.set_Matrix(matrix);
				gr.set_Matrix(m);
			}
			//dVertex.clear();
		}
	}
	//dVertex.clear();
	for (auto& thr : threads) {
		thr.join();
	}
	//for (int i = 0; i < threads.size(); i++) {
	//	threads_graph[i].detach();
	//}
}

void check_k_connected(std::vector<Graph> arr) {
	int numVector = 0;
	if (!arr.empty()) {
		int k = arr[0].get_k();
		for (auto v : arr) {
			std::cout << "Matrix number: " << numVector << std::endl;
			Matrix m = v.get_Matrix();
			m.PrintMatrix();
			int numEdge = m.currentNumEdge();
			std::cout << "numEdge " << numEdge << std::endl;
			std::cout << "CONGRATULATE!!!! It is k-connected graph; k = " << k << std::endl;
			std::cout << std::endl;
			numVector++;
		}
	}
	std::cout << "END!" << std::endl;
}

void check() {
	int numNode = 100;
	int numEdge = 0;
	int k = 7;
	int numLine = numNode;
	int numColumn = numNode;
	Graph graph(numNode, k);
	unsigned int start_time = clock();
	bool res = false;
	res = graph.algorithmEven();
	unsigned int end_time = clock();
	unsigned int search_time = end_time - start_time;
	std::cout << "Time work = " << search_time / CLOCKS_PER_SEC << std::endl;
	res = graph.checkMinGraph();
	unsigned int end_time2 = clock();
	unsigned int search_time2 = end_time2 - end_time;
	std::cout << "Time work = " << search_time2 / CLOCKS_PER_SEC << std::endl;
	res = graph.checkContractionMinmality();
	unsigned int end_time3 = clock();
	unsigned int search_time3 = end_time3 - end_time2;
	std::cout << "Time work = " << search_time3 / CLOCKS_PER_SEC << std::endl;
}

void printToFile(std::vector<Graph> graphs, std::string nameFile) {
	std::fstream f;
	f.open(nameFile, std::ios::out);
	if (f)
	{
		std::string str;
		int number = 1;
		for (auto f_graph : graphs)
		{
			int numLine = f_graph.get_numLine();
			int numColumn = f_graph.get_numColumn();
			Matrix m = f_graph.get_Matrix();
			str = "N ";
			str += std::to_string(number);
			f << str << std::endl;
			str = "";
			for (int i = 0; i < numLine; i++) {
				str = "T ";
				str += std::to_string(i);
				f << str << std::endl;
				str = "";
			}
			for (int i = 0; i < numLine; i++) {
				for (int j = i + 1; j < numColumn; j++) {
					if (m.getElem(i, j) == 1) {
						str = "E ";
						str += std::to_string(i);
						str += " ";
						str += std::to_string(j);
						str += " ";
						str += std::to_string(m.getElem(i, j));
						f << str << std::endl;
						str = "";
					}
				}
			}
			number++;
			f << std::endl;
		}
		f.close();
	}
}

int main(void) {
	int numNode = 7;
	int numEdge = 0;
	int k = 4;
	int numLine = numNode;
	int numColumn = numNode;
	Graph graph(numNode, k);
	unsigned int start_time = clock();
	std::cout << "Start algorithm " << std::endl;
	std::cout << std::endl;
	bool res;
	int num = 0;
	std::vector<Graph> arr_graph;
	
	/*Matrix matr(numLine, numColumn);
	matr.setElem(0, 3, 1);
	matr.setElem(0, 4, 1);
	matr.setElem(0, 5, 1);
	matr.setElem(0, 6, 1);
	matr.setElem(1, 2, 1);
	matr.setElem(1, 4, 1);
	matr.setElem(1, 5, 1);
	matr.setElem(1, 6, 1);
	matr.setElem(2, 1, 1);
	matr.setElem(2, 3, 1);
	matr.setElem(2, 5, 1);
	matr.setElem(2, 6, 1);
	matr.setElem(3, 0, 1);
	matr.setElem(3, 2, 1);
	matr.setElem(3, 4, 1);
	matr.setElem(3, 6, 1);
	matr.setElem(4, 0, 1);
	matr.setElem(4, 1, 1);
	matr.setElem(4, 3, 1);
	matr.setElem(4, 5, 1);
	matr.setElem(5, 0, 1);
	matr.setElem(5, 1, 1);
	matr.setElem(5, 2, 1);
	matr.setElem(5, 4, 1);
	matr.setElem(6, 0, 1);
	matr.setElem(6, 1, 1);
	matr.setElem(6, 2, 1);
	matr.setElem(6, 3, 1);
	graph.set_Matrix(matr);*/

	std::cout << std::endl;
	res = graph.algorithmEven();
	if (res) {
		arr_graph.push_back(graph);
		array_graph.push_back(graph);
	}
	//unsigned int end_time = clock();
	//unsigned int search_time = end_time - start_time;
	//std::cout << "Time work = " << search_time / CLOCKS_PER_SEC << std::endl;
	//..fun(graph, &num, k, &arr_graph);
	//Thread_Graph th;
	//th.thread_graph(graph, k);
	//arr_graph = th.get_array();
	thread_graph(graph, k);
	int numVector = 0;

	//check_k_connected(arr_graph);	
	arr_graph = array_graph;
	std::cout << std::endl;

	//std::cout << "The end!" << std::endl;

	bool r = false;
	std::string str = "";
	std::vector<Graph> min_graph;
	for (auto v : arr_graph) {
		r = false;
		r = v.checkMinGraph();
		//std::cout << std::endl;
		if (r == true) {
			v.set_minimality(r);
			min_graph.push_back(v);
		}
	}

	std::vector<Graph> finish_graph;
	for (auto m_graph : min_graph) {
		r = false;
		//m_graph.get_Matrix().PrintMatrix();
		r = m_graph.checkContractionMinmality();
		//std::cout << std::endl;
		if (r == true) {
			m_graph.set_minContraction(r);
			finish_graph.push_back(m_graph);
		}
	}
	for (auto f_graph : finish_graph) {
		std::cout << "Matrix number: " << numVector << std::endl;
		Matrix m = f_graph.get_Matrix();
		m.PrintMatrix();
		int numEdge = m.currentNumEdge();
		std::cout << "numEdge " << numEdge << std::endl;
		std::cout << "CONGRATULATE!!!! It is k-connected graph; k = " << k << std::endl;
		std::cout << "Result of check min graph = true" << std::endl;
		std::cout << "Result of check min contraction graph = true" << std::endl;
		std::cout << std::endl;
		numVector++;
	}
	std::cout << "Num " << num + 1 << std::endl;
	std::cout << std::endl;
	unsigned int end_time = clock();
	unsigned int search_time = end_time - start_time;
	std::cout << "Time work = " << search_time / CLOCKS_PER_SEC << std::endl;
	printToFile(finish_graph, "out1.txt");
	arr_graph.clear();
	min_graph.clear();
	finish_graph.clear();

	//check();
	return 0;
}