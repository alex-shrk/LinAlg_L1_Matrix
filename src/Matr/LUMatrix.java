package Matr;

public class LUMatrix extends Matrix {
    private double determinant = 0.;

    public LUMatrix(int sizeSquareMatr) {
        super(sizeSquareMatr);
    }

    public double getElemL(int row, int col) {
        if (row > col) {
            return matrix[row][col];
        } else if (row == col)
            return 1.0;
        return 0.0;
    }

    public double getElemU(int row, int col) {
        if (row <= col) {
            return matrix[row][col];
        }
        return 0.0;
    }

    public void setElemL(int row, int col, double val) {
        if (row > col)
            matrix[row][col] = val;


    }

    public void setElemU(int row, int col, double val) {
        if (row <= col) {
            matrix[row][col] = val;
        }
    }

    public void calculateDeterminant() {

        double result = 1.0;
        for (int i = 0; i < m; i++) {
            if (matrix[i][i] == 0) {
                result = matrix[i][i];
                break;
            }
            result *= matrix[i][i];
        }
        this.determinant = result;


    }

    public double getDeterminant() {
        if (this.determinant==0.)
            calculateDeterminant();
        return this.determinant;
    }


    public Matrix getLMatrix(){
        Matrix l=new Matrix(this.m);
        for(int i=0;i<m;i++){
            for (int j=0;j<n;j++)
                l.matrix[i][j]=getElemL(i,j);
        }
        return l;
    }
    public Matrix getUMatrix(){
        Matrix u=new Matrix(this.m);
        for(int i=0;i<m;i++){
            for (int j=0;j<n;j++)
                u.matrix[i][j]=getElemU(i,j);
        }
        return u;
    }


    public Matrix getInverseMatrix(){//todo проверка u and optimization
        Matrix inverseMatr=getEyeMatrix(m);


        for (int i=m-1;i>=0;i--){
            for (int j=m-1;j>=0;j--){
                double val_x=0.;
                if (i==j){
                    val_x=1.;
                    for (int k=j+1;k<m;k++){
                        val_x -=getElemU(j,k)*inverseMatr.matrix[k][j];
                    }
                    inverseMatr.matrix[i][j]=val_x/getElemU(i,j);
                }

                else{
                    if (i<j){
                        val_x=0.;
                        for (int k=i+1;k<m;k++){
                            val_x -=getElemU(i,k)*inverseMatr.matrix[k][j];
                        }
                        inverseMatr.matrix[i][j]=val_x/getElemU(i,i);
                    }
                    else {
                        val_x=0.;
                        for (int k=j+1;k<m;k++){
                            val_x -=inverseMatr.matrix[i][k]*getElemL(k,j);
                        }
                        inverseMatr.matrix[i][j]=val_x;

                    }
                }



            }
        }
        return inverseMatr;
    }



}
