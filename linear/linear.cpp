#include <iostream>
#include <cmath>
using namespace std;
#define eps1 pow(10, -3)
#define eps2 pow(10, -6)
// размерность матрицы от 4 до 6

double* gauss(double** matrixx, double* col, int n) // метод Гаусса
{
    double* result = new double[n];
    int i, j, l;
    double boo;
    for (l = 0; l < n; l++) // прямой ход
    {
        if (matrixx[l][l] == 0)
        {
            for (int t = l+1; t < n; t++)
            {
                if (matrixx[t][l] != 0) swap(matrixx[t], matrixx[l]);
            }
        }
        for (i = l + 1; i < n; i++)
        {
            boo = matrixx[i][l];
            for (j = l; j < n; j++)
            {
                matrixx[i][j] -= (matrixx[l][j] / matrixx[l][l]) * boo;
            }
            col[i] -= (col[l] / matrixx[l][l]) * boo;
        }
    }
    result[n-1] = col[n-1]/matrixx[n-1][n-1];
    
    for (i = n - 2; i >= 0; i--) // обратный
    {
        result[i] = col[i];
        for (j = n - 1; j > i; j--)
        {
            result[i] -= matrixx[i][j] * result[j];
        }
        result[i] /= matrixx[i][i];
    }
    return result;
}

double* multiply(double** matrixx, double* col, int n)
{
    double* result = new double[n];
    for (int i = 0; i < n; i++)
    {
        result[i] = 0;
        for (int j = 0; j < n; j++)
        {
            result[i] += matrixx[i][j] * col[j];
        }
    }
    return result;
}

double** multiply(double** matrix, int n)
{
    double** result = new double* [n];
    int i, j, k;
    for (k = 0; k < n; k++)
    {
       result[k] = new double[n];
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            result[i][j] = 0;
            for (k = 0; k < n; k++)
            {
                result[i][j] += matrix[i][k]*matrix[j][k];
            }
        }
    }
    return result;
}

double** root(double** matrixx, int n) // нахождение матрицы S (верхнетреугольной)
{
    double** triangular = new double* [n];
    int i, j, k;
    for (k = 0; k < n; k++)
    {
        triangular[k] = new double[n];
    }
    triangular[0][0] = sqrt(matrixx[0][0]);
    for (i = 1; i < n; i++)
    {
        triangular[0][i] = matrixx[0][i] / triangular[0][0];
    }
    double sum;
    for (i = 1; i < n; i++)
    {
        sum = 0;
        for (j = 0; j < i; j++)
        {
            sum += pow(triangular[j][i], 2);
        }
        triangular[i][i] = sqrt(matrixx[i][i] - sum);
        for (j = 0; j < i; j++) 
        {
            triangular[i][j] = 0;
        }
        for (j = i + 1; j < n; j++)
        {
            sum = 0;
            for (k = 0; k < i; k++)
            {
                sum += triangular[k][i] * triangular[k][j];
            }
            triangular[i][j] = (matrixx[i][j] - sum) / triangular[i][i];
        }
    }
    return triangular;
} 

double* solve(double** matrix, double* col, int n, string k) // решаем системы с треугольными матрицами
{
    double* result = new double[n];
    int i, j;
    if (k == "down") // если нижнетреугольная
    {
        result[0] = col[0] / matrix[0][0];
        for (i = 1; i < n; i++)
        {
            result[i] = col[i];
            for (j = 0; j < i; j++)
            {
                result[i] -= matrix[i][j] * result[j];
            }
            result[i] /= matrix[i][i];
        }
        return result;
    }
    if (k == "up") // если верхнетреугольная
    {
        result[n - 1] = col[n - 1] / matrix[n - 1][n - 1];

        for (i = n - 2; i >= 0; i--) 
        {
            result[i] = col[i];
            for (j = n - 1; j > i; j--)
            {
                result[i] -= matrix[i][j] * result[j];
            }
            result[i] /= matrix[i][i];
        }
        return result;
    }

}

void trans(double** matrix, int k)
{
    double** temp = new double* [k];
    for (int i = 0; i < k; i++)
    {
        temp[i] = new double[k];
    }
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
        {
            temp[i][j] = matrix[j][i];
        }
    }
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
        {
            matrix[i][j] = temp[i][j];
        }
    }
    for (int i = 0; i < k; i++)
    {
        delete[] temp[i];
    }
    delete[] temp;
}

int main()
{
    setlocale(0, "");
    int method;
    cout << "Введите номер требуемого метода: 1 - Гаусс, 2 - квадратный корень, 3 - простой итерационный:" << endl;
    cin >> method;
    for (int dim = 3; dim <= 6; dim++) // начинаем с СЛАУ размерности три (которая для проверки), 
                                       // потом плохо обусловленные от размерности 4 до 6
    {
        if ((dim > 3)&&(method == 2)) return 0;
        int i, j, k;
        double** matrix = new double* [dim]; // исходная матрица, потом будет меняться в методе Гаусса 
        double** matrix0 = new double* [dim]; // матрица-дублер, для проверок
        double** Bmatrix = new double* [dim]; // для метода простой итерации
        double* column = new double[dim]; // столбец свободных членов
        double* res; // результат, иначе вектор х
        const double N = 19; // номер варианта
        for (k = 0; k < dim; k++)
        {
            matrix[k] = new double[dim];
            matrix0[k] = new double[dim];
            Bmatrix[k] = new double[dim];
        }
        if (dim == 3)
        {
            for (k = 0; k < dim; k++)
            {
                for (int u = 0; u < dim; u++)
                {
                    if (u == k)
                    {
                        matrix[k][u] = N + pow(2, u + 1);
                        matrix0[k][u] = matrix[k][u];
                    }
                    else
                    {
                        matrix[k][u] = 1;
                        matrix0[k][u] = matrix[k][u];
                    }
                }
            }
            for (k = 0; k < dim; k++)
            {
                column[k] = N + 2 * (k + 2);
            }
        }
        else
        {
            double eps;
            if ((method == 1) || (method == 2)) eps = eps1;
            else eps = eps2;
            for (i = 0; i < dim; i++)
            {
                for (j = 0; j < dim; j++)
                {
                    if (i == j)
                    {
                        matrix[i][j] = 1 + eps * N;
                        matrix0[i][j] = matrix[i][j];
                    }
                    else 
                    {
                        if (i > j)
                        {
                            matrix[i][j] = eps * N;
                            matrix0[i][j] = matrix[i][j];
                        }
                        else
                        {
                            matrix[i][j] = -1 - eps * N;
                            matrix0[i][j] = matrix[i][j];
                        }
                    }
                }
            }
            for (i = 0; i < dim - 1; i++)
            {
                column[i] = -1;
            }
            column[dim - 1] = 1;
        }
        if (method == 1) // гаусса
        {
            double* res = gauss(matrix, column, dim);
            cout << "Результат:";
            for (k = 0; k < dim; k++)
            {
                cout << res[k] << endl;
            }
            cout << "Проверка:";
            res = multiply(matrix0, res, dim);
            for (k = 0; k < dim; k++)
            {
                cout << res[k] << endl;
            }
        }
        if (method == 2) // квадратного корня
        {
            double** choles = root(matrix0, dim);
            for (int i = 0; i < dim; i++) // выводим матрицу S
            {
                for (int j = 0; j < dim; j++)
                {
                    cout << choles[i][j] << " ";
                }
                cout << endl;
            }
            res = solve(choles, column, dim, "up"); // решаем систему Sy = b
            for (int i = 0; i < dim; i++)
            {
                cout << res[i] << endl;
            }
            trans(choles, dim); // получим нижнетреугольную матрицу
            res = solve(choles, res, dim, "down"); // решаем систему St*x = y
            cout << "Результат:" << endl;
            for (int i = 0; i < dim; i++)
            {
                cout << res[i] << endl;
            }
            cout << endl << "Проверка:" << endl;
            res = multiply(matrix0, res, dim);
            for (k = 0; k < dim; k++)
            {
                cout << res[k] << endl;
            }
        }
        if (method == 3)
        {
            double normB = 0;
            double maxnorm = 0;
            for (i = 0; i < dim; i++) // строим матрицу B, а также вычислим ее норму (inf)
            {
                if (matrix0[i][i] == 0)
                {
                    for (k = i + 1; k < dim; k++)
                    {
                        if (matrix0[k][i] != 0) swap(matrix0[i], matrix[k]);
                    }
                }
                normB = 0;
                for (j = 0; j < dim; j++)
                {
                    Bmatrix[i][j] = -matrix0[i][j] / matrix0[i][i];
                    if (i == j) Bmatrix[i][j] = 0;
                    cout << Bmatrix[i][j] << " ";
                    normB += abs(Bmatrix[i][j]);
                }
                cout << endl;
                if (normB > maxnorm) maxnorm = normB;
                column[i] /= matrix0[i][i];
            }
            cout << "Норма матрицы B: " << maxnorm << endl;
            double* x1 = new double[dim];
            for (i = 0; i < dim; i++)
            {
                x1[i] = column[i];
            }
            double* x2 = new double[dim];
            double norm;
            if (maxnorm >= 1) cout << "Итерационный процесс неприменим!" << endl;
            else
            {
                do
                {
                    norm = 0;
                    x2 = multiply(Bmatrix, x1, dim);
                    for (i = 0; i < dim; i++)
                    {
                        x2[i] += column[i];
                        norm += abs(x2[i] - x1[i]);
                        x1[i] = x2[i];
                    }
                } while (maxnorm * norm / (1 - maxnorm) > eps2);
                cout << "Результат:" << endl;
                for (i = 0; i < dim; i++)
                {
                    cout << x2[i] << endl;
                }
                cout << "Проверка:" << endl;
                res = multiply(matrix0, x2, dim);
                for (i = 0; i < dim; i++)
                {
                    cout << res[i] << endl;
                }
            }
        }
        for (i = 0; i < dim; i++)
        {
            delete[] matrix[i];
            delete[] matrix0[i];
            delete[] Bmatrix[i];
        }
        delete[] matrix;
        delete[] matrix0;
        delete[] Bmatrix;
        delete[] column;
    }
    return 0;
}
 
