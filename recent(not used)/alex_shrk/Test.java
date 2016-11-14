import Exceptions.MatrixSizeError;
import Matr.Matrix;

public class Test {
    public static void main(String[] args) throws MatrixSizeError {
        Matrix a=new Matrix(new double[][]{{5.,1.,-3.},{3.,0.,-2.},{-4.,-1.,1.}});
        //Matrix a=new Matrix(new double[][]{{3./5.,0.,4./5.},{0.,1.,0.},{-4./5.,0.,3./5.}});
        //Matrix a=new Matrix(new double[][]{{3.,-1.,0.},{0.,0.,0.},{4.,7.,0.}});
        /*QRMatrix qr=Matrix.decompositionQRRotare(a);
        Matrix q=qr.getQ();
        q.outputToTxt("q_rot");
        Matrix r=qr.getR();

        r.outputToTxt("r_rot");
        q.multiply(r).outputToTxt("a_");*/
        a.eigenvaluesMatrixQR(0.001).outputToTxt("eigen");
        a.eigenvaluesMatrixQR_old(0.001).outputToTxt("eig_old");

        //Matrix a=new Matrix(3);
        //a.generateMatrix();
        /*QRMatrix P=Matrix.hessenbergMatrix(a);
        Matrix R=P.getQ().transpose();
        Matrix a_=R.multiply(a).multiply(R.transpose());
        a_.outputToTxt("H_");*/
        //Matrix.decompositionQRRotare(a);


    }
}
