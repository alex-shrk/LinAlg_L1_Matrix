package Matr;

public interface IMatrix {
    int getM();
    int getN();
    int getCountColumns();
    int getCountRows();
    double getElement(int m, int n);

    void setElement(int m, int n, double el) throws ArrayIndexOutOfBoundsException;

    void outputToTxt(String path);




}
