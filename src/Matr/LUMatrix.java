package Matr;

public class LUMatrix extends Matrix {
    private double determinant = 0.;

    public LUMatrix(int sizeSquareMatr) {
        super(sizeSquareMatr);
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
                double value_x=0.;
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



}
