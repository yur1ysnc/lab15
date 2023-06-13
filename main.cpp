#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
class Matrix {
private:
    int rows, cols;
    double** data;
public:
    Matrix() : rows(0), cols(0), data(nullptr) {}
	Matrix(int rows_t,int cols_t) {
		rows = rows_t;
		cols = cols_t;
		data = new double*[rows];
		for (int i = 0; i < rows; ++i) {
            data[i] = new double[cols];
      
		}
	}
	Matrix(const std::string& filename) {
        std::ifstream file(filename);
        if (file.is_open()) {
            file >> rows >> cols;
           
            data = new double*[rows];
				for(int i = 0; i < rows; i++) {
				data[i] = new double[cols];
					for(int j = 0; j < cols; j++)
						file >> data[i][j];
			}
			file.close();
        } 
		else {
            throw std::runtime_error("Failed to open file");
        }
    }
    Matrix(int rows_t, int cols_t, double** data2) {
    	rows = rows_t;
		cols = cols_t;
        data = new double*[rows];
        for (int i = 0; i < rows; ++i) {
            data[i] = new double[cols];
            for (int j = 0; j < cols; ++j) {
                data[i][j] = data2[i][j];
            }
        }
    }
    Matrix(const Matrix& other)  {
    	rows = other.rows;
		cols = other.cols;
        data = new double*[rows];
        for (int i = 0; i < rows; ++i) {
            data[i] = new double[cols];
            for (int j = 0; j < cols; ++j) {
                data[i][j] = other.data[i][j];
            }
        }
    }
    ~Matrix() {
        if (rows > 0){
            for (int i = 0; i < rows; ++i) {
                delete[] data[i];
            }
            delete[] data;
	}
    }
    Matrix& operator=(const Matrix& other) {
        if (&other == this) {
            return *this;
        }
        if (rows > 0){
            for (int i = 0; i < rows; ++i) {
                delete[] data[i];
            }
            delete[] data;
	}
        rows = other.rows;
        cols = other.cols;
        data = new double*[rows];
        for (int i = 0; i < rows; ++i) {
            data[i] = new double[cols];
            for (int j = 0; j < cols; ++j) {
                data[i][j] = other.data[i][j];
            }
        }
        return *this;
    }
    double* operator[](int row) {
        return data[row];
    }
    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw std::invalid_argument("Matrices have incompatible dimensions");
        }
		double** temp_ar;
		temp_ar = new double*[rows];
	for (int i = 0; i < rows;++i){
	    temp_ar[i] = new double[other.cols];
	    for (int j = 0; j < other.cols;++j){
	        temp_ar[i][j] = 0;		
	    }
	}
        Matrix result(rows, other.cols,temp_ar);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
	        double sum = double();
                for (int k = 0; k < cols; ++k) {
                    sum += data[i][k] * other.data[k][j];
                }
                result[i][j] = sum;
            }
        }
        return result;
	for (int i = 0; i < 2;++i){delete[] temp_ar[i];}
	delete[] temp_ar;
    }
    Matrix operator+(const Matrix& other) const {
        if (cols != other.cols || rows != other.cols) {
            throw std::invalid_argument("Matrices have incompatible dimensions");
        }
		double** temp_ar;
		temp_ar = new double*[rows];
		for (int i = 0; i < rows;++i){
			temp_ar[i] = new double[cols];
			for (int j = 0; j < cols;++j){
			temp_ar[i][j] = 0;	
			}
		}
		Matrix result(rows,cols,temp_ar);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j){
			result[i][j] = data[i][j] + other.Getdata()[i][j];	
			}
		}
	    	for (int i = 0; i < 2;++i){delete[] temp_ar[i];}
		delete[] temp_ar;
		return result;
	}
	Matrix get_transp() const {
		if (rows!=cols) {throw std::invalid_argument("Matrix has incompatible dimension");}
		double** temp_ar;
		temp_ar = new double*[rows];
		for (int i = 0; i < rows;++i){
			temp_ar[i] = new double[cols];
			for (int j = 0; j < cols;++j){
			temp_ar[i][j] = data[j][i];	
			}
		}
		Matrix result(rows,cols,temp_ar);
		return result;
	}
    Matrix operator-(const Matrix& other) const {
        if (cols != other.cols || rows != other.cols) {
            throw std::invalid_argument("Matrices have incompatible dimensions");
        }
		double** temp_ar;
		temp_ar = new double*[rows];
		for (int i = 0; i < rows;++i){
			temp_ar[i] = new double[cols];
			for (int j = 0; j < cols;++j){
			temp_ar[i][j] = 0;	
			}
		}
		Matrix result(rows,cols,temp_ar);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j){
			result[i][j] = data[i][j] - other.Getdata()[i][j];	
			}
		}
		for (int i = 0; i < 2;++i){delete[] temp_ar[i];}
		delete[] temp_ar;
		return result;
		
    }
    Matrix operator*(double scalar) const {
	double** temp_ar;
	temp_ar = new double*[rows];
	for (int i = 0; i < rows;++i){
		temp_ar[i] = new double[cols];
		for (int j = 0; j < cols;++j){
			temp_ar[i][j] = 0;
		}
	}
        Matrix result(rows, cols,temp_ar);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result[i][j] = data[i][j] * scalar;
            }
        }
	for (int i = 0; i < 2;++i){delete[] temp_ar[i];}
	delete[] temp_ar;
        return result;
    }
	Matrix operator/(double scalar) const {
		double** temp_ar;
		temp_ar = new double*[rows];
		for (int i = 0; i < rows;++i){
			temp_ar[i] = new double[cols];
			for (int j = 0; j < cols;++j){
				temp_ar[i][j] = 0;
			}
		}
		Matrix result(rows, cols,temp_ar);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result[i][j] = data[i][j] / scalar;
			}
		}
		for (int i = 0; i < 2;++i){delete[] temp_ar[i];}
		delete[] temp_ar;
		return result;
		}
    friend Matrix operator*(const double& scalar, const Matrix& matrix) {return matrix * scalar;}
    int Getrows() const{
        return rows;
    }
    int Getcols() const{
	return cols;
    }
	double GetDet() const {
		
		if (rows != cols) {throw std::runtime_error("Matrix should be square, error!");}
		else {
			std::vector<int> ex_cols = {};
			return p_GetDet(0, rows - 1, rows - 1, ex_cols);
		}
	}
	double p_GetDet(int row_first, int row_last, int col_last, std::vector<int>& ex_cols) const {
		if (row_first == row_last) {
			for(int j = 0; j <= col_last; j++) {
				if(std::count(ex_cols.begin(),ex_cols.end(),j) == 0)
					return data[row_first][j];
			}
		}
		double result  = 0;
		int col_parity = 0;
		for(int j = 0; j <= col_last; j++) {
			if(std::count(ex_cols.begin(),ex_cols.end(),j) == 0) {
				ex_cols.push_back(j);
				result += ((col_parity++ % 2 == 0) ? 1 : -1) * this->data[row_first][j] * p_GetDet(row_first + 1, row_last, col_last, ex_cols);
				ex_cols.pop_back();
			}
		}
		return result;
	}
	Matrix operator!() const {
		if(rows != cols) {
            std::cerr << "Matrix should be square, error!" << std::endl;
            exit(1);
        }
        return get_Inv();
    }
	Matrix get_Inv() const {
		double d = GetDet();
		if (d==0){
			std::cerr << "Impossible to create inverse matrix because of zero determinate";
			exit(1);
		}
		Matrix   tr(get_transp());
		Matrix   result(tr);
		for(int i = 0; i < result.rows; i++) {
			for(int j = 0; j < result.rows; j++) {
				Matrix   lower(rows-1,cols-1);
				int l_i = 0;
				for(int ii = 0; ii < result.rows; ii++) {
					int l_j = 0;
					bool flag = false;
					for(int jj = 0; jj < result.rows; jj++) {
						if(ii != i && jj != j) {
							flag = true;
							lower.data[l_i, l_j++] = tr.data[ii, jj];
						}
					}
					if(flag)
						l_i++;
				}
				result.data[i][j] =  lower.GetDet() * (((i + j) % 2 == 0) ? 1 : -1);
			}
		}
		return result / d;

		}
    double** Getdata() const{return data;}
    void writedoubleoFile(const std::string& filename) const {
        std::ofstream file(filename);
        if (file.is_open()) {
            file << rows << " " << cols << std::endl;
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    file << data[i][j] << " ";
                }
                file << std::endl;
            }
            file.close();
        } else {
            throw std::runtime_error("Failed to open file");
        }
    }
	static Matrix ZeroCr(int rows,int cols){
		int r = rows;
		int c = cols;
		double** temp_ar_l;
		temp_ar_l = new double*[r];
		for (int i = 0; i < r;++i){
			temp_ar_l[i] = new double[c];
			for (int j = 0; j < c;++j){
				temp_ar_l[i][j] = 0;
			}
		}
		Matrix result_l(r, c,temp_ar_l);
		for (int i = 0; i < 2;++i){delete[] temp_ar_l[i];}
		delete[] temp_ar_l;
		return result_l;
		
	}
	static Matrix UnitCr(int rows,int cols){
		int r = rows;
		int c = cols;
		double** temp_ar_l;
		temp_ar_l = new double*[r];
		for (int i = 0; i < r;++i){
			temp_ar_l[i] = new double[c];
			for (int j = 0; j < c;++j){
				if (i==j) {temp_ar_l[i][j] = 1;}
				else {temp_ar_l[i][j] = 0;}
			}
		}
		Matrix result_l(r, c,temp_ar_l);
		for (int i = 0; i < 2;++i){delete[] temp_ar_l[i];}
		delete[] temp_ar_l;
		return result_l;
	}
	bool operator==(const Matrix& other) {
		int c = 0;
		if ((rows != other.rows) || (cols != other.cols) ){c+=1;}
		else {
			for (int i = 0; i < rows; ++i ) {
				for (int j = 0; j < cols;++j){
					if (data[i][j] != other.data[i][j]){
						c+=1;
					}
			    }	
			}
		}
		if (c==0) {return true;}
		return false;
	}
    friend std::ostream& operator<<(std::ostream &os,const Matrix& matr){
        for (int i = 0; i < matr.Getrows();++i){
			for (int j = 0; j < matr.Getcols();++j){
				os << matr.Getdata()[i][j] << " ";
			}
			os << "\n";
		}
		return os;
    }
};
int main() {
    double **data2;
    double **data3;
    int k = 2;
    int l = 4;
    data2 = new double*[2];
    for (int i = 0; i < 2; ++i){
        data2[i] = new double[2];
        for (int j = 0; j < 2; ++j){
            data2[i][j] = k;
            k++;
        }
    }
    data3 = new double*[2];
    for (int i = 0; i < 2; ++i){
        data3[i] = new double[2];
        for (int j = 0; j < 2; ++j){
            data3[i][j] = l;
            l++;
        }
    }
    Matrix matr2(2,2,data2);
    Matrix matr3(2,2,data3);
    std::cout <<"matr2 = " << "\n" <<  matr2 << std::endl;// '<<' demonstration 
	std::cout << "matr2 determinat  = " << matr2.GetDet() << std::endl;
	std::cout << "matr2 transp = \n " << matr2.get_transp() << std::endl;
	std::cout << "!matr2  = \n " << (!matr2) << std::endl;
    std::cout <<"matr3 = " << "\n" <<  matr3 << std::endl;// '<<' demonstration 
    Matrix e = matr2;
    std::cout <<"e = " << "\n" <<  e << std::endl;    // ' = ' demonstration
    std::cout <<"matr2 == matr3 (0 - no, 1 - yes) -> " <<  (matr2 == matr3) << std::endl; 
    std::cout <<"matr2 == e (0 - no, 1 - yes) -> " <<  (matr2 == e) << std::endl; // '==' demonstration
    std::cout <<"matr3*matr2 = " << "\n" <<  matr3*matr2 << std::endl;// '*' demonstration 
    std::cout <<"matr3+matr2 = " << "\n" <<  matr3+matr2 << std::endl;// '+' demonstration
    std::cout <<"matr3-matr2 = " << "\n" <<  matr3-matr2 << std::endl; // '-' demonstration 
    std::cout <<"matr3 * 5 =  " << "\n" <<  matr3*5 << std::endl;   // '* scalar' demonstration
    std::cout <<"Zero 3x3 matrix = \n" <<  Matrix::ZeroCr(3,3) << std::endl; // ' zero matrix ' demonstration
    std::cout <<"doublehe unit 3x3 matrix = \n" <<  Matrix::UnitCr(3,3) << std::endl; // ' unit matrix ' demonstration  
}