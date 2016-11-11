import Exceptions.MatrixSizeError;
import Matr.Matrix;
import Matr.QRMatrix;

public class Test {
    public static void main(String[] args) throws MatrixSizeError {
        Matrix a=new Matrix(new double[][]{{3.,-1.,5.},{0.,0.,5.},{4.,7.,5.}});
        QRMatrix qr=a.decompositionQR();
        qr.getQ().outputToTxt("q");
        qr.getR().outputToTxt("r");

    }
}
