package Matr;

public class QRMatrix{
    protected Matrix q;
    protected Matrix r;
    private int m;//size QR Matrix

    public QRMatrix(int sizeQRMatrix) {
        this.m = sizeQRMatrix;
        this.q=new Matrix(m);
        this.r=new Matrix(m);
    }
    protected void setValuesOfMatrix(double val){
        for (int i=0;i<m;i++)
            for(int j=0;j<m;j++) {
                q.matrix[i][j] = val;
                r.matrix[i][j] = val;
            }
    }

    public Matrix getQ() {
        return q;
    }

    public Matrix getR() {
        return r;
    }

    public int getM() {
        return m;
    }
}
