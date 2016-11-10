import Exceptions.MatrixSizeError;
import Matr.LUMatrix;
import Matr.Matrix;
import Matr.QRMatrix;

import java.util.Scanner;

public class Lab1 {
    public static void main(String[] args) {
        try {
            int MAX_COUNT_OF_ITERATIONS=1000000;
            Scanner sc=new Scanner(System.in);
            System.out.println("Введите размерность матрицы: ");
            int size_matr = sc.nextInt();


            System.out.println("Введите степень точности для нахождения собственных значений с помощью qr - разложения");
            double precision_qr = Matrix.getPrecision(sc.nextInt());



            for (int exc=1;exc<=2;exc++){
                Matrix matrix = null;
                String typeMatrix = null;
                if (exc==1){
                    typeMatrix="goodCondMatrix";
                    System.out.println("Выполнение операций для хорошообусловленной матрицы");

                    for (int i=0;i<MAX_COUNT_OF_ITERATIONS;i++){
                        Matrix m =new Matrix(size_matr);
                        //m.setEps(degree_eps);
                        m.generateMatrix();
                        if (m.conditionMatrix(m.decompLU())<1.0e8) {
                            matrix = m;
                            break;
                        }
                    }
                    assert matrix != null;
                    matrix.outputToTxt(typeMatrix);
                }
                else
                    if (exc==2){
                        typeMatrix="badCondMatrix";
                        System.out.println("Выполнение операций для плохообусловленной матрицы");
                        matrix=Matrix.getBadCondMatrix(size_matr);
                        matrix.outputToTxt(typeMatrix);
                    }

                System.out.println("1.");
                Matrix x = Matrix.generateVector(matrix.getM());
                Matrix b=matrix.multiply(x);//Ax=b
                System.out.println(" b:");
                b.outputToScreen();
                b.outputToTxt("b_"+typeMatrix);

                System.out.println("2.");
                LUMatrix luMatrix = matrix.decompLU();
                Matrix solveX = matrix.solveSystem(b,luMatrix);
                System.out.println("Solve X:");
                solveX.transpose().outputToScreen();
                solveX.transpose().outputToTxt("solveX_"+typeMatrix);

                System.out.println("3.");
                System.out.println("    Determinant A = "+luMatrix.getDeterminant());
                System.out.println("    Condition A = "+matrix.conditionMatrix(luMatrix));

                System.out.println("5.");

                double delta = x.subtraction(solveX).sphericalNorm()/x.sphericalNorm();
                System.out.println("    Относительная погрешность решения СЛАУ: "+delta);


                System.out.println("6.");


                QRMatrix QR = matrix.decompQR();
                Matrix AAt = matrix.multiply(matrix.transpose());
                System.out.println("    Собственные значения матрицы AAt:");
                Matrix eigenvaluesMatrixQR = AAt.eigenvaluesMatrixQR(precision_qr);
                eigenvaluesMatrixQR.outputToScreen();
                eigenvaluesMatrixQR.outputToTxt("eigenvaluesMatrixQR_"+typeMatrix);

                double eigenvalueMax=Double.MIN_VALUE;
                double eigenvalueMin=Double.MAX_VALUE;
                for (int i=0;i<size_matr;i++){
                    eigenvalueMax=Math.max(eigenvalueMax,Math.abs(eigenvaluesMatrixQR.getElement(i)));
                    eigenvalueMin=Math.min(eigenvalueMin,Math.abs(eigenvaluesMatrixQR.getElement(i)));
                }
                double eigenvalCheck=Math.sqrt(eigenvalueMax/eigenvalueMin);
                System.out.println("sqrt(eigenvalue Max/eigenvalue Min)= "+eigenvalCheck);
                System.out.println();
                System.out.println("__________________________________________________");
                System.out.println();



            }









        } catch (MatrixSizeError matrixSizeError) {
            matrixSizeError.printStackTrace();
        }


    }
}
