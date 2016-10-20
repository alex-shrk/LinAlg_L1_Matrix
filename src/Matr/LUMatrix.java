package Matr;

import Exceptions.MatrixSizeError;

public class LUMatrix extends Matrix {
    public LUMatrix(int sizeSquareMatr) {
        super(sizeSquareMatr);
    }
    public double getElemL(int line,int col){
        if (line>col) {
            return matrix[line][col];
        }
        else if(line==col)
            return 1.0;
        return 0.0;
    }
    public double getElemU(int line,int col){
        if (line<=col){
            return matrix[line][col];
        }
        return 0.0;
    }

    public void setElemL(int line,int col,double val){
        if (line>col)
            matrix[line][col]=val;
        else if (line==col){}

    }
    public void setElemU(int line,int col,double val){
        if (line<=col){
            matrix[line][col]=val;
        }
    }
    public double determinant() {
        double result = 1.0;
        for (int i = 0; i < m; i++) {
            result *= matrix[i][i];
        }
        return result;
    }
    public LUMatrix(Matrix A) throws MatrixSizeError {
        int m=A.getM();
        int n=A.getN();
        if (m!=n){
            throw new MatrixSizeError("Матрица не квадратная");
        }
        LUMatrix lu=new LUMatrix(m);
        lu.initializeValuesOfMatr(0.0);

        for (int i=0;i<m;i++){
            for (int j=i;j<n;j++){
                double val_u=A.getElement(i,j);
                for (int k=0;k<i;k++) {
                    val_u -= lu.getElemL(i, k) * lu.getElemU(k, j);
                }
                lu.setElemU(i,j,val_u);

                double val_l=A.getElement(j,i);
                for (int k=0;k<i;k++){
                    val_l -= lu.getElemL(j,k) * lu.getElemU(k,i);
                }
                val_l/=lu.getElemU(i,i);
                lu.setElemL(j,i,val_l);

            }
        }


    }






}
