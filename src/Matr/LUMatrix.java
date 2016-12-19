package Matr;

import Exceptions.MatrixSizeError;

public class LUMatrix extends Matrix {
    private double determinant = 0.;

    public LUMatrix(int sizeSquareMatr) {
        super(sizeSquareMatr);
    }
    public LUMatrix(IMatrix matr) throws MatrixSizeError {
        m=matr.getCountRows();
        int n=matr.getCountColumns();

        if (m!=n){
            throw new MatrixSizeError("Матрица не является квадратной"+m+"!="+n);
        }

        this.matrix = new double[m][n];

        for (int i=0;i<m;i++){
            for (int j=i;j<n;j++){
                double value_u=matr.getElement(i,j);//a[i][index]

                for (int k=0;k<i;k++) {
                    value_u -= getElemL(i,k) * getElemU(k,i);

                }
                setElemU(i,j,value_u);

                if (j==i) continue;

                double value_l=matr.getElement(j,i);
                for (int k=0;k<i;k++){
                    value_l -= getElemL(j,k) *getElemU(k,i) ;

                }

                //if (Math.abs(getElemU(i,i))>=eps) {
                    value_l /= getElemU(i,i);
                //}
                //else{
                //    System.out.println("Элемент u "+i+" "+i+" меньше eps");
                //    break;
                //}
                setElemL(j,i,value_l);


            }
        }


    }

    public LUMatrix(IMatrix matr,boolean fillindMode,int startValueFilling) throws MatrixSizeError {
        m=matr.getCountRows();
        int n=matr.getCountColumns();

        if (m!=n){
            throw new MatrixSizeError("Матрица не является квадратной"+m+"!="+n);
        }

        this.matrix = new double[m][m];

        int counter_of_filling = 0;

        for (int i=0;i<m;i++){
            for (int j=i;j<m;j++){
                double value_u=matr.getElement(i,j);//a[i][index]

                for (int k=0;k<i;k++) {
                    value_u -= getElemL(i,k) * getElemU(k,i);

                }
                setElemU(i,j,value_u);

                if (matr.getElement(i,j)==0.0d && value_u!=0.0d)
                    counter_of_filling++;

                if (j==i) continue;

                double value_l=matr.getElement(j,i);
                for (int k=0;k<i;k++){
                    value_l -= getElemL(j,k) *getElemU(k,i) ;

                }

                if (Math.abs(getElemU(i,i))>=eps) {
                    value_l /= getElemU(i,i);
                }
                else{
                    System.out.println("Элемент u "+i+" "+i+" меньше eps");
                    break;
                }
                setElemL(j,i,value_l);
                if (matr.getElement(j,i)==0.0d && value_l!=0.0d)
                    counter_of_filling++;

            }
        }
        if (fillindMode)
            System.out.println("Заполнение для метода LU = "+(counter_of_filling-startValueFilling));
    }


    public double getElemL(int row, int column) {
        if (row == column) {
            return 1.0;
        }
        else {
            if (row > column) {
                return matrix[row][column];
            }
        }
        return 0.0;
    }

    public double getElemU(int row, int column) {
        if (row <= column) {
            return matrix[row][column];
        }
        return 0.0;
    }

    public void setElemL(int row, int column, double value) {
        if (row > column)
            matrix[row][column] = value;


    }

    public void setElemU(int row, int column, double value) {
        if (row <= column) {
            matrix[row][column] = value;
        }
    }

    public void calculateDeterminant() {

        double result = 1.0;

        for (int i = 0; i < m; i++) {
            if (matrix[i][i] == 0) {
                result = 0.;
                break;
            }
            result *= matrix[i][i];
        }
        this.determinant = result;

    }

    public double getDeterminant() {
        if (this.determinant==0.) {
            calculateDeterminant();
        }
        return this.determinant;
    }


    public Matrix getLMatrix(){
        Matrix l=new Matrix(this.m);
        for(int i=0;i<m;i++){
            for (int j=0;j<n;j++) {
                l.matrix[i][j] = getElemL(i, j);
            }
        }
        return l;
    }
    public Matrix getUMatrix(){
        Matrix u=new Matrix(this.m);
        for(int i=0;i<m;i++){
            for (int j=0;j<n;j++) {
                u.matrix[i][j] = getElemU(i, j);
            }
        }
        return u;
    }



    public Matrix getInverseMatrix(){
        Matrix inverseMatr=getEyeMatrix(m);


        for (int i=m-1;i>=0;i--){

            for (int j=m-1;j>=0;j--){
                double value_x=1.;
                if (i==j){
                    value_x=1.;
                    for (int k=j+1;k<m;k++){
                        value_x -=getElemU(j,k)*inverseMatr.matrix[k][j];
                    }
                    inverseMatr.matrix[i][j]=value_x/getElemU(i,j);

                }
                else{
                    if (i<j){
                        value_x=0.;
                        for (int k=i+1;k<m;k++){
                            value_x -=getElemU(i,k)*inverseMatr.matrix[k][j];
                        }
                        inverseMatr.matrix[i][j]=value_x/getElemU(i,i);

                    }
                    else {
                        value_x=0.;
                        for (int k=j+1;k<m;k++){
                            value_x -=inverseMatr.matrix[i][k]*getElemL(k,j);
                        }
                        inverseMatr.matrix[i][j]=value_x;

                    }
                }

            }
        }
        return inverseMatr;

    }

    public Matrix solveSystem(IMatrix b) throws MatrixSizeError {
        if (m != b.getCountRows()) {
            throw new RuntimeException("Число строк матрицы не совпадает с числом строк матрицы b");
        }

        Matrix x = new Matrix(m, 1);
        Matrix y = new Matrix(m, 1);

        for (int i = 0; i < m; i++) {
            double value = b.getElement(i, 0);
            for (int k = 0; k < i; k++) {
                value -= getElemL(i, k) * y.getElement(k);
            }
            y.setElement(i, value);
        }

        for (int i = m - 1; i >= 0; i--) {
            double value = y.getElement(i);
            for (int k = i + 1; k < m; k++) {
                value -= getElemU(i, k) * x.getElement(k);
            }
            x.setElement(i, value / getElemU(i, i));
        }

        return x;


    }



}
