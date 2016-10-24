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
    protected double eps=0.0001;//todo сделать нормальное представление
    private LUMatrix luMatrix;

    private static int counterLU=0;




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


    /*public static Matrix inverse(Matrix matr) throws MatrixSizeError, MatrixNotExist {
        int m=matr.getCountRows();
        int n=matr.getCountColumns();
        Matrix inversedMatrix=new Matrix(m,n);

        double det=determinant();
        if (det==0)
            throw new MatrixNotExist("Обратная матрица не существует");
        Matrix invMatr=new Matrix(m,n);
        for (int i=m-1;i>=0;i--){

        }


    }*/
/*

    */
/** метод обрежет матрицу до размеров int order*//*

    public static Matrix cropMatrix(int order,Matrix A){
        Matrix B=new Matrix(order);
        for (int i=0;i<order;i++){
            System.arraycopy(A.matrix[i],0,B.matrix[i],0,order);
        }
        return B;
    }
*/

    /*public static double getMajorMinor(int order,Matrix A) throws Exception {
        //order-порядок главного минора матрицы А
        if (order>A.matrix.length)
            throw new Exception("порядок минора больше порядка матрицы");

        return determinant(cropMatrix(order, A));
    }*/

    /*public static boolean checkingPositivityOfMajorMinors(Matrix A) throws Exception {
        for (int i=1;i<A.matrix.length;i++) {
            if (getMajorMinor(i, A) <= 0) {
                return false;
            }
        }
        return true;

    }*/

    /*public static double determinant(Matrix A) throws MatrixSizeError {
        if (A.getM()!=A.getN())
            throw new MatrixSizeError("Матрица не квадратная");
        double det=1.0;
        List<Matrix> lu=A.decompLU(A);
        for (int i=0;i<lu.get(1).getCountColumns();i++)
            det*=lu.get(1).matrix[i][i];
        return det;
    }*/







    //Сферическая векторная норма

    /*public double norma2() {
        Matrix AA;
        if (m <= n) {
            AA = this.multiply(transpose());
        } else {
            AA = transpose().multiply(this);
        }

        Double max = 0.0;
        for (int k = 0; k < AA.m; k++) {
            max = Math.max(max, AA.data[k][k]);
        }
        return Math.sqrt(max);
    }*/



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

    public int getCountOfOperationsLUDecomp(){
        return counterLU;
    }
/*

    public static LUMatrix decompLULead(Matrix matr,double eps) throws MatrixSizeError, CloneNotSupportedException {
        int m_matr=matr.getCountRows();
        int n_matr=matr.getCountColumns();
        double[][] matr_copy= Arrays.copyOf(matr.matrix,m_matr);

        int[] q = new int[m_matr];
        for (int i=0;i<m_matr;i++) {//initialize value of q vector
            q[i] = i;
            matr_copy[i] = Arrays.copyOf(matr.matrix[i], m_matr);
        }

        if (m_matr!=n_matr){
            throw new MatrixSizeError("Матрица не квадратная");
        }
        LUMatrix luMatrix=new LUMatrix(m_matr);
        luMatrix.setValuesOfMatrix(0.);

        //transposition columns
        for (int k=0;k<m_matr;k++){
            int m=0;
            for (int i=0;i<m_matr;i++)
                m= luMatrix.getElemU(k,m)<luMatrix.getElemU(k,i) ? i : m;
            if (m!=0){
                double tmp_m;
                int tmp_q;
                for (int i=0;i<m_matr;i++){
                    tmp_m=matr_copy[i][k];
                    matr_copy[i][k]=matr_copy[i][m];
                    matr_copy[i][m]=tmp_m;
                }
                tmp_q=q[k];
                q[k]=q[m];
                q[m]=tmp_q;

            }


            for (int i=0;i<m_matr;i++) {
                for (int j = i; j < n_matr; j++) {
                    double val_u = matr_copy[i][j];//a[i][j]
                    for (int k1 = 0; k1 < i; k1++) {
                        val_u -= luMatrix.matrix[i][k1] * luMatrix.matrix[k1][j];
                    }
                    luMatrix.matrix[i][j] = val_u;

                    if (j == i) continue;
                    double val_l = matr_copy[j][i];
                    for (int k1 = 0; k1 < i; k1++) {
                        val_l -= luMatrix.matrix[j][k1] * luMatrix.matrix[k1][i];
                    }
                    //if (Math.abs(luMatrix.matrix[i][i]) >= eps)
                        val_l /= luMatrix.matrix[i][i];
                    */
/*else {
                        System.out.println("Элемент u" + i + " " + i + "меньше eps");
                        break;
                    }*//*



                    luMatrix.matrix[j][i] = val_l;

                }
            }

        }



        return luMatrix;

    }
    public LUMatrix decompLULead() throws MatrixSizeError, CloneNotSupportedException {
        this.luMatrixLead= decompLULead(this,eps);
        return this.luMatrixLead;
    }
*/


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


    /*public Matrix inverse() {
        double[][] result = new double[m][m];

        for (int p = m-1; p >= 0; p--) {
            Double val = 1.0;
            for (int k = p+1; k < m; k++) {
                val -= getElemU(p, k) * result[k][p];
            }
            val /= getElemU(p, p);
            result[p][p] = val;

            for (int i = p-1; i >= 0; i--) {
                val = 0.0;
                for (int k = i+1; k < m; k++) {
                    val -= getElemU(i, k) * result[k][p];
                }
                val /= getElemU(i, i);
                result[i][p] = val;
            }

            for (int j = p-1; j >= 0; j--) {
                val = 0.0;
                for (int k = j+1; k < m; k++) {
                    val -= result[p][k] * getElemL(k, j);
                }
                result[p][j] = val;
            }
        }

        return new Matrix(result);
    }*/





    /*public void inputFromFile(String path) throws FileNotFoundException {
            System.out.println("Inputing from file...");
            File file = new File(path);
            if (!file.exists()){
                throw new FileNotFoundException(file.getName());
            }
            try{
               // BufferedReader in = new BufferedReader(new FileReader(file.getAbsoluteFile()));
                Scanner sc=new Scanner(file);
                sc.useDelimiter(" ");
                try{
                    while (sc.hasNextDouble()){
                        String[] value=line.split(" ");
                        for (int i = 0; i < n; i++) {
                            for (int j = 0; j < n; j++) {
                                matrix[i][j] = Double.parseDouble(value[i]);
                            }
                        }
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
                finally {
                    in.close();
                }
            }
            catch (IOException e){
                throw new RuntimeException(e);
            }
            System.out.println("Complete");

        }*/

    public  void outputToTxt(String path){
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
