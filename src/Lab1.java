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
            System.out.println("Введите размерность хорошо обусловленной матрицы: ");
            int size_good_matr = sc.nextInt();
            System.out.println("Введите размерность плохо обусловленной матрицы: ");
            int size_bad_matr = sc.nextInt();
            double goodCond=1.e8;


            System.out.println("Введите точность для нахождения собственных значений с помощью qr - разложения");
            double precision_qr = sc.nextDouble();



            for (int var=1;var<=2;var++){
                Matrix matrix = null;
                String typeMatrix = null;
                if (var==1){
                    typeMatrix="goodCondMatrix";
                    System.out.println("Выполнение операций для хорошо обусловленной матрицы\n");

                    for (int i=0;i<MAX_COUNT_OF_ITERATIONS;i++){
                        Matrix m =new Matrix(size_good_matr);
                        m.generateMatrix();
                        if (m.conditionMatrix()<goodCond) {
                            matrix = m;
                            break;
                        }
                    }
                    assert matrix != null;
                    matrix.outputToTxt(typeMatrix);
                }
                else
                    if (var==2){
                        typeMatrix="badCondMatrix";
                        System.out.println("Выполнение операций для плохо обусловленной матрицы\n");
                        matrix=Matrix.getBadCondMatrix(size_bad_matr);
                        matrix.outputToTxt(typeMatrix);
                    }
                
                System.out.println("Задание 1.\n");

                System.out.println("Создание СЛАУ Ax=b");
                Matrix x = Matrix.generateVector(matrix.getM());
                Matrix b=matrix.multiply(x);//Ax=b
                System.out.println("Вектор b:");
                b.outputToTxt("b_"+typeMatrix);

                System.out.println("Задание 2.\n");
                System.out.println("Решение СЛАУ с помощью LU-разложения");
                LUMatrix luMatrix = matrix.decompositionLU();
                Matrix solveX = matrix.solveSystem(b,luMatrix);
                System.out.println("Решение X:");
                solveX.outputToTxt("solveX_"+typeMatrix);

                System.out.println("Задание 3.\n");
                System.out.println("Вычисление определителя и числа обусловленности для матрицы СЛАУ");
                System.out.println("    Determinant A = "+luMatrix.getDeterminant());
                System.out.println("    Condition A = "+matrix.conditionMatrix());

                System.out.println("\nЗадание 5.\n");
                System.out.println("Вычисление относительной погрешности решения СЛАУ");
                double d1=x.subtraction(solveX).euclideanNorm();
                double d2=x.euclideanNorm();
                double delta = d1/d2;
                System.out.println("    delta=: "+delta);


                System.out.println("Задание 6.\n");
                System.out.println("Нахождение собственных значений матрицы AAt с помощью QR-алгоритма");


                QRMatrix QR = matrix.decompositionQR();
                Matrix AAt = matrix.multiply(matrix.transpose());
                System.out.println("    Собственные значения матрицы AAt:");
                Matrix eigenvaluesMatrixQR = AAt.eigenvaluesMatrixQR(precision_qr);
                eigenvaluesMatrixQR.outputToTxt("eigenvaluesMatrixQR_"+typeMatrix);
                double eigenvaluesCheck=Matrix.checkEigenvaluesMatrix(eigenvaluesMatrixQR);
                System.out.println("sqrt(max(eigenvalue) /min(eigenvalue))= "+eigenvaluesCheck);
                System.out.println();
                System.out.println("__________________________________________________");
                System.out.println();


            }
        } catch (MatrixSizeError matrixSizeError) {
            matrixSizeError.printStackTrace();
        }


    }
}
