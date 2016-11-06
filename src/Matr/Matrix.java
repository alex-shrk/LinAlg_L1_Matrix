package Matr;

import Exceptions.MatrixSizeError;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.GregorianCalendar;
import java.util.Random;

public class Matrix {
    protected double[][] matrix;
    //m-count of rows
    //n-count of columns
    protected int m,n;
    protected double rangeMin=0.0;
    protected double rangeMax=10.0;

    protected double eps=0.0001;
    private LUMatrix luMatrix;

    private static int counterLU=0;

    public void setEps(int degree){
        this.eps=1/Math.pow(10,degree);
    }




    public Matrix() {
    }




    public int getM() {
        return m;
    }

    public int getN() {
        return n;
    }

    public void setRange(double rangeMin,double rangeMax){
        this.rangeMin=rangeMin;
        this.rangeMax=rangeMax;
    }
    public Matrix(int m, int n) {
        this.m=m;
        this.n=n;
        this.matrix = new double[m][n];
    }
    public Matrix(int sizeSquareMatr){
        this(sizeSquareMatr,sizeSquareMatr);
    }

    public Matrix(double[][] values){
        this(values.length,values[0].length);
        for (int i=0;i<m;i++){
            System.arraycopy(values[i],0,this.matrix[i],0,n);
        }
    }
    public static Matrix getEyeMatrix(int sizeMatr){
        Matrix matr=new Matrix(sizeMatr);
        for (int i=0;i<sizeMatr;i++)
            matr.matrix[i][i]=1.;
        return matr;
    }
    protected void setValuesOfMatrix(double val){
        for (int i=0;i<m;i++)
            for(int j=0;j<n;j++)
                matrix[i][j]=val;
    }

    public int getCountColumns() {
        return n;
    }

    public int getCountRows() {
        return m;
    }

    public double getElement(int m,int n) throws ArrayIndexOutOfBoundsException {
        try{
            return matrix[m][n];
        }
        catch (Exception e) {
            throw new ArrayIndexOutOfBoundsException("Element" + m + ' ' + n + "not found");
        }
    }
    public double getElement(int i) throws ArrayIndexOutOfBoundsException {
        try{
            if (m==1)
                return matrix[0][i];
            else
            if (n==1)
                return matrix[i][0];
            else throw new MatrixSizeError("Невозможно получить элемент матрицы. Только для векторов");
        }
        catch (Exception e) {
            throw new ArrayIndexOutOfBoundsException("Element" + m + ' ' + n + "not found");
        }
    }
    public void setElement(int m,int n,double el) throws ArrayIndexOutOfBoundsException {
        try{
            matrix[m][n]=el;
        }
        catch (Exception e) {
            throw new ArrayIndexOutOfBoundsException("Element" + m + ' ' + n + "not found");
        }
    }
    public void setElement(int i, double el) throws ArrayIndexOutOfBoundsException {
        try{
            if (m==1)
                matrix[0][i]=el;
            else
            if (n==1)
                matrix[i][0]=el;
            else throw new MatrixSizeError("Невозможно присвоить элемент матрице. Только для векторов");
        }
        catch (Exception e) {
            throw new ArrayIndexOutOfBoundsException("Element" + m + ' ' + n + "not found");
        }
    }
    public Matrix multiply(Matrix matr) throws MatrixSizeError {
        return multiply(this,matr);
    }

    public static Matrix multiply(Matrix A, Matrix B) throws MatrixSizeError {
        int m_a=A.getCountRows();
        int n_a=A.getCountColumns();

        int m_b=B.getCountRows();
        int n_b=B.getCountColumns();

        if (n_a!=m_b)
            throw new MatrixSizeError("Произведение матриц невозможно: число столбцов 1-й матрицы!=числу строк 2-й матрицы");

        Matrix c=new Matrix(m_a,n_b);
        for(int i=0;i<c.m;i++){
            for(int j=0;j<c.n;j++){
                c.matrix[i][j]=0.0;
                for (int k=0;k<n_a;k++)
                    c.matrix[i][j]+=A.matrix[i][k]*B.matrix[k][j];
            }
        }
        return c;
    }

    public static Matrix multiply(Matrix a, double val){
        int m=a.getCountRows();
        int n=a.getCountColumns();

        for(int i=0;i<m;i++){
            for(int j=0;j<n;j++){
                a.matrix[i][j]*=val;
            }
        }
        return a;
    }
    public Matrix multiply(double val){
        return multiply(this,val);
    }
    public static Matrix addition(Matrix A, Matrix B) throws MatrixSizeError {
        int n_a=A.getCountColumns();
        int n_b=B.getCountColumns();

        int m_a=A.getCountRows();
        int m_b=B.getCountRows();

        if (m_a!=m_b && n_a!=n_b)
            throw new MatrixSizeError();
        double el_a,el_b;

        Matrix C=new Matrix(m_a,n_a);

        for(int i=0;i<m_a;i++){
            for(int j=0;j<n_a;j++){
                el_a=A.getElement(i,j);
                el_b=B.getElement(i,j);
                C.matrix[i][j]=el_a+el_b;
            }
        }
        return C;
    }
    public Matrix addition(Matrix B) throws MatrixSizeError {
        return addition(this,B);
    }
    public static Matrix subtraction(Matrix A, Matrix B) throws MatrixSizeError {
        int n_a=A.getCountColumns();
        int n_b=B.getCountColumns();

        int m_a=A.getCountRows();
        int m_b=B.getCountRows();

        if (m_a!=m_b && n_a!=n_b)
            throw new MatrixSizeError();
        double el_a,el_b;

        Matrix C=new Matrix(m_a,n_a);

        for(int i=0;i<m_a;i++){
            for(int j=0;j<n_a;j++){
                el_a=A.getElement(i,j);
                el_b=B.getElement(i,j);
                C.matrix[i][j]=el_a-el_b;
            }
        }
        return C;
    }
    public Matrix subtraction(Matrix B) throws MatrixSizeError{
        return subtraction(this,B);
    }

    /**
     * @return
     * ( A 0 )
     * ( 0 B )
     */
    public static Matrix combineMatrix(Matrix A,Matrix B){
        int m=A.getM()+B.getM();
        int n=A.getN()+B.getN();

        Matrix C = new Matrix(m,n);
        C.setValuesOfMatrix(0.);

        int m_a=A.getM();

        for (int i=0;i<m_a;i++){
            System.arraycopy(A.matrix[i],0,C.matrix[i],0,A.getN());
        }
        for (int i=m_a,j=0;i<m;i++,j++){
            System.arraycopy(B.matrix[j],0,C.matrix[i],A.getN(),B.getN());
        }
        return C;
    }


    public double getRandomDoubleValue(){
        Random random = new Random();
        return rangeMin+(rangeMax-rangeMin)*random.nextDouble();
    }
    public boolean isCorrectRange(double rangeMin,double rangeMax){
        return !Double.valueOf(rangeMax - rangeMin).isInfinite();
    }
    public boolean isCorrectRange(){
        return isCorrectRange(rangeMin,rangeMax);
    }



    public void generateMatrix() {
        System.out.println("Generating matrix. Size="+m+"x"+n);
        if (isCorrectRange()) {

            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    matrix[i][j] = getRandomDoubleValue();
                }
            }
            System.out.println("Complete");
        }
        else throw new ArithmeticException("rangeMax-rangeMin=infinite");
    }

    public Matrix generateVector(){
        return generateVector(this.m);
    }
    public static Matrix generateVector(int m){
        int n=1;
        Matrix x=new Matrix(m,n);
            for (int i=0;i<m;i++){
                x.setElement(i,i+1.);
            }

        return x;
    }

    public static Matrix transpose(Matrix matr){
        int m=matr.getCountRows();
        int n=matr.getCountColumns();
        Matrix transposedMatrix=new Matrix(n,m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                transposedMatrix.matrix[j][i] = matr.matrix[i][j];
            }
        }
        return transposedMatrix;

    }
    public Matrix transpose(){
        return transpose(this);
    }


    public static Matrix getPartOfMatrix(int startIndexM,int endIndexM,int startIndexN,int endIndexN,Matrix A){
        assert (startIndexM<endIndexM)&&(startIndexN<endIndexN);

        int m=endIndexM-startIndexM+1;
        int n=endIndexN-startIndexN+1;

        Matrix B=new Matrix(m,n);
        for (int i=startIndexM;i<=endIndexM;i++){
            System.arraycopy(A.matrix[i],startIndexN,B.matrix[i-startIndexM],0,n);
        }
        return B;
    }
    public Matrix getPartOfMatrix(int startIndexM,int endIndexM,int startIndexN,int endIndexN){
        return getPartOfMatrix(startIndexM,endIndexM,startIndexN,endIndexN,this);
    }


    /**
     * euclidean norm
    */
    public static double sphericalNorm(Matrix a) {
        int m=a.getM();
        int n=a.getN();
        double norm=0.0;
        for (int i=0;i<m;i++){
            for (int j=0;j<n;j++){
                norm+=Math.pow(Math.abs(a.matrix[i][j]),2);
            }
        }
        return Math.sqrt(norm);
    }
    public double sphericalNorm(){
        return sphericalNorm(this);
    }

    public static double conditionMatrix(Matrix a,LUMatrix lu){
        double sphericalNormA=sphericalNorm(a);
        double sphericalNormInvA=sphericalNorm(lu.getInverseMatrix());
        return sphericalNormA*sphericalNormInvA;
    }
    public double conditionMatrix(){
        return conditionMatrix(this,this.luMatrix);
    }


    public LUMatrix decompLU(Matrix matr) throws MatrixSizeError {
        this.counterLU=0;
        int m=matr.getCountRows();
        int n=matr.getCountColumns();

        if (m!=n){
            throw new MatrixSizeError("Матрица не квадратная");
        }
        LUMatrix luMatrix=new LUMatrix(m);
        luMatrix.setValuesOfMatrix(0.);


        for (int i=0;i<m;i++){
            for (int j=i;j<n;j++){
                double val_u=matr.matrix[i][j];//a[i][j]
                for (int k=0;k<i;k++) {
                    val_u -= luMatrix.matrix[i][k] * luMatrix.matrix[k][j];
                    counterLU+=2;
                }
                luMatrix.matrix[i][j]=val_u;

                if (j==i) continue;
                double val_l=matrix[j][i];
                for (int k=0;k<i;k++){
                    val_l -= luMatrix.matrix[j][k] * luMatrix.matrix[k][i];
                    counterLU+=2;
                }
                if (Math.abs(luMatrix.matrix[i][i])>=eps) {
                    val_l /= luMatrix.matrix[i][i];
                    counterLU++;
                }
                else{
                    System.out.println("Элемент u"+i+" "+i+"меньше eps");
                    break;
                }


                luMatrix.matrix[j][i]=val_l;

            }
        }


        return luMatrix;

    }

    public static QRMatrix decompQR(Matrix a) throws MatrixSizeError{
        int m=a.getCountRows();
        int n=a.getCountColumns();

        if (m!=n){
            throw new MatrixSizeError("Матрица не квадратная");
        }
        QRMatrix qrMatrix=new QRMatrix(m);
        qrMatrix.setValuesOfMatrix(0.);
        Matrix p;
        Matrix h;
        Matrix h1 = null;

        for (int k=0;k<m-1;k++){
            p=a.getPartOfMatrix(k,m-1,k,k);
            double el=p.getElement(0);
            double norm=p.sphericalNorm();
            p.setElement(0, el + (el >=0 ? norm : -norm));
            h = p.multiply(p.transpose()).multiply(2.0/(p.transpose().multiply(p).getElement(0)));
            h=Matrix.getEyeMatrix(m-k).subtraction(h);
            h=Matrix.combineMatrix(getEyeMatrix(k),h);
            h1 = (h1!=null) ? h.multiply(h1) : h;

        }

        assert h1 != null;
        qrMatrix.r=h1.multiply(a);
        qrMatrix.q=h1.transpose();
        return qrMatrix;
    }
    public QRMatrix decompQR() throws MatrixSizeError{
        return decompQR(this);
    }






    public LUMatrix decompLU() throws MatrixSizeError{
        this.luMatrix= decompLU(this);
        return this.luMatrix;
    }


    public static Matrix solveSystem(Matrix matrix,LUMatrix luMatrix, double eps) throws MatrixSizeError {
        int m=matrix.getM();
        Matrix y=new Matrix(m,1);
        y.setValuesOfMatrix(0.0);


        Matrix x=generateVector(m);
        Matrix b=generateVector(m);

        //Matrix b=new Matrix(new double[][]{{1.}, {2.}, {3.},{4.}});

        double u=0.0;
        double sum=0.0;
        y.matrix[0][0]=b.matrix[0][0];

        for (int i=1;i<=m-1;i++) {
            sum=0.;
            for (int k = 0; k < i; k++) {
                sum += luMatrix.getElemL(i, k) * y.matrix[k][0];
                counterLU+=2;
            }
            y.matrix[i][0] = b.matrix[i][0] - sum;
            counterLU++;

        }
        for (int i=m-1;i>=0;i--){
            sum=0.;
            for(int k=i+1;k<m;k++){
                sum+=luMatrix.getElemU(i,k)*x.matrix[k][0];
                counterLU+=2;
            }
            u=luMatrix.getElemU(i,i);
            if (Math.abs(u)>=eps) {
                x.matrix[i][0] = (y.matrix[i][0] - sum) / u;
                counterLU+=2;
            }

            else{
                System.out.println("Элемент u"+i+" "+i+"меньше eps");
                break;
            }

        }
        System.out.println("Число операций при решении методом LU: "+counterLU);
        counterLU=0;

        return x;


    }
    public Matrix solveSystem(LUMatrix luMatrix) throws MatrixSizeError {
        return solveSystem(this,luMatrix,eps);
    }




    public void outputToTxt(String path){
        System.out.println("Outputing to file...");
        String default_path="matrix"+ new GregorianCalendar().getTime().toString();
        if (path.isEmpty())
            path=default_path;
        File file = new File(path+".txt");
        try{
            if (!file.exists())
                file.createNewFile();

            PrintWriter out=new PrintWriter(file.getAbsoluteFile());
            try{

                for (int i = 0; i < m; i++) {
                    for (int j = 0; j < n; j++) {
                        out.print(matrix[i][j]);
                        out.print('\t');
                    }
                    out.print('\n');
                }
            }
            finally {
                out.close();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("Complete");
    }


}
