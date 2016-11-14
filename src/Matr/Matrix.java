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
    protected double rangeMin=-10.0;
    protected double rangeMax=10.0;

    protected double eps=1.e-5;
    private static int counterSumSub=0;
    private static int counterMultDivSum=0;



    public double getEps() {
        return eps;
    }
    public void setEps(int eps){
        this.eps=eps;
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



    /** return vector (1,2,...m)*/
    public static Matrix generateVector(int m){
        int n=1;
        Matrix x=new Matrix(m,n);
        for (int i=0;i<m;i++){
            x.setElement(i,i+1.);
        }

        return x;
    }

    public Matrix generateVector(){
        return generateVector(this.m);
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
    public static Matrix getBadCondMatrix(int size_matr){
        Matrix badCondMatr=new Matrix(size_matr);
        for (int i=0;i<size_matr;i++) {
            for (int j = i; j < size_matr; j++) {
                if (j == i) {
                    badCondMatr.setElement(j, j, (0.01 / ((size_matr - j + 1) * (j + 1))));
                    continue;
                }

                badCondMatr.setElement(i, j, 0.);//upper matrix
                badCondMatr.setElement(j, i, i * (size_matr - j));//lower matrix
            }
        }
        return badCondMatr;
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

    /** Only for vector*/
    public double getElement(int i) throws ArrayIndexOutOfBoundsException {
        try{
            if (m==1) {
                return matrix[0][i];
            }

            else {

                if (n == 1) {
                    return matrix[i][0];
                }
                else
                    throw new MatrixSizeError("Невозможно получить элемент матрицы. Операция применима только к векторам");
            }

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
            if (m==1){
                matrix[0][i]=el;
            }

            else{

                if (n==1) {
                    matrix[i][0] = el;
                }

                else throw new MatrixSizeError("Невозможно получить элемент матрицы. Операция применима только к векторам");

            }
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
            throw new MatrixSizeError("Невозможно выполнить операцию умножения матриц: число столбцов 1-й матрицы!=числу строк 2-й матрицы");

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
        int m_a=A.getCountRows();
        int m_b=B.getCountRows();

        int n_a=A.getCountColumns();
        int n_b=B.getCountColumns();


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
                norm+=Math.abs(a.matrix[i][j]*a.matrix[i][j]);
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
    public double conditionMatrix(LUMatrix LUMatrix){
        return conditionMatrix(this,LUMatrix);
    }


    public static double checkEigenvaluesMatrix(Matrix eigenvaluesMatrix){
        double eigenvalueMax=Double.MIN_VALUE;
        double eigenvalueMin=Double.MAX_VALUE;
        for (int i=0;i<eigenvaluesMatrix.m;i++){
            eigenvalueMax=Math.max(eigenvalueMax,Math.abs(eigenvaluesMatrix.getElement(i)));
            eigenvalueMin=Math.min(eigenvalueMin,Math.abs(eigenvaluesMatrix.getElement(i)));
        }
        System.out.println("Max "+eigenvalueMax);
        System.out.println("Min "+eigenvalueMin);
        return Math.sqrt(eigenvalueMax/eigenvalueMin);
    }


    public LUMatrix decompositionLU(Matrix matr) throws MatrixSizeError {
        counterSumSub=0;
        counterMultDivSum=0;
        int m=matr.getCountRows();
        int n=matr.getCountColumns();

        if (m!=n){
            throw new MatrixSizeError("Матрица не является квадратной"+m+"!="+n);
        }

        LUMatrix LUMatrix=new LUMatrix(m);
        LUMatrix.setValuesOfMatrix(0.);


        for (int i=0;i<m;i++){
            for (int j=i;j<n;j++){
                double value_u=matr.matrix[i][j];//a[i][j]

                for (int k=0;k<i;k++) {
                    value_u -= LUMatrix.matrix[i][k] * LUMatrix.matrix[k][j];
                    counterSumSub++;
                    counterMultDivSum++;
                }
                LUMatrix.matrix[i][j]=value_u;

                if (j==i) continue;

                double value_l=matrix[j][i];
                for (int k=0;k<i;k++){
                    value_l -= LUMatrix.matrix[k][i] * LUMatrix.matrix[j][k];
                    counterSumSub++;
                    counterMultDivSum++;
                }
                if (Math.abs(LUMatrix.matrix[i][i])>=eps) {
                    value_l /= LUMatrix.matrix[i][i];
                    counterMultDivSum++;
                }
                else{
                    System.out.println("Элемент u "+i+" "+i+" меньше eps");
                    break;
                }
                LUMatrix.matrix[j][i]=value_l;

            }
        }


        return LUMatrix;

    }

    public static QRMatrix decompositionQRRotare(Matrix a) throws MatrixSizeError{
        int m=a.getCountRows();
        int n=a.getCountColumns();
        double eps=1.e-2;

        if (m!=n){
            throw new MatrixSizeError("Матрица не квадратная");
        }
        QRMatrix qrMatrix=new QRMatrix(m);
        Matrix p;
        //Matrix H;
        Matrix Q=null;
        Matrix R=a;


        Matrix G;
        int i_ind;
        int j_ind;


        for (int k=0;k<m;k++){
            p = R.getPartOfMatrix(0, m-1, k, k);
            for (int j=k;j<p.m;j++) {
                double v = p.getElement(j);
                if (Math.abs(v) <= eps)//элемент уже нулевой с определенной точностью
                    continue;
                //i=m-1 j=j
                if (p.m-1>j) {
                    i_ind = j;
                    j_ind = p.m-1;
                }
                else{
                    i_ind = p.m-1;
                    j_ind = j;
                }
                double a_i=p.getElement(i_ind);
                double a_j=p.getElement(j_ind);
                double norm = Math.sqrt((a_i*a_i)+(a_j*a_j));
                double cos=a_i/norm;
                double sin=a_j/norm;

                G=Matrix.getEyeMatrix(m);
                G.setElement(i_ind,i_ind,cos);
                G.setElement(j_ind,j_ind,cos);
                G.setElement(i_ind,j_ind,sin);
                G.setElement(j_ind,i_ind,-sin);


                if (Q != null) {
                    Q = G.multiply(Q);
                } else {
                    Q = G;
                }
                R = Q.multiply(a);
            }
        }
        assert Q != null;
        qrMatrix.q=Q.transpose();
        qrMatrix.r=R;
        return qrMatrix;
    }

    public static QRMatrix decompositionQR(Matrix a) throws MatrixSizeError{
        int m=a.getCountRows();
        int n=a.getCountColumns();

        if (m!=n){
            throw new MatrixSizeError("Матрица не квадратная");
        }
        QRMatrix qrMatrix=new QRMatrix(m);
        Matrix p;
        Matrix H;
        Matrix Q=null;
        Matrix R=a;


        for (int k=0;k<m-1;k++){

            p=R.getPartOfMatrix(k,m-1,k,k);
            double v=p.getElement(0);//v*e_1
            double norm=p.sphericalNorm();
            if (v<0)//sign(v)
                norm=-norm;
            p.setElement(0, v + norm);
            Matrix ptp=p.transpose().multiply(p);
            H = p.multiply(p.transpose()).multiply(2./ptp.getElement(0));

            H=Matrix.getEyeMatrix(m-k).subtraction(H);
            H=Matrix.combineMatrix(getEyeMatrix(k),H);

            if (Q!=null) {
                Q = H.multiply(Q);
            }
            else {
                Q = H;
            }
            R=Q.multiply(a);
        }
        assert Q != null;
        qrMatrix.q=Q.transpose();
        qrMatrix.r=R;
        return qrMatrix;
    }
    public QRMatrix decompositionQR() throws MatrixSizeError{
        return decompositionQR(this);
    }

    /*public static Matrix hessenbergMatrix(Matrix A) throws MatrixSizeError {
        int m=A.getCountRows();
        int n=A.getCountColumns();
        double eps=1.e-8;

        if (m!=n){
            throw new MatrixSizeError("Матрица не квадратная");
        }
//        if (m<=n+1)
//            throw new MatrixSizeError("Приведение к матрице Хессенберга невозможно.  "+m+">"+n+"+1")

        Matrix p;
        Matrix H;
        Matrix Q=null;
        Matrix R=A;

        //применяем метод отражений Хаусхолдера и зануляем нижние элементы матрицы лежащие ниже на 2 диагонали
        for (int k=0;k<m-2;k++){

            p=R.getPartOfMatrix(k+1,m-1,k,k);
            double v=p.getElement(0);//v*e_1
            double norm=p.sphericalNorm();
            if (v<0)//sign(v)
                norm=-norm;
            p.setElement(0, v + norm);
            Matrix ptp=p.transpose().multiply(p);
            H = p.multiply(p.transpose()).multiply(2./ptp.getElement(0));

            H=Matrix.getEyeMatrix(m-(k+1)).subtraction(H);
            H=Matrix.combineMatrix(getEyeMatrix(k+1),H);

            if (Q!=null) {
                Q = H.multiply(Q);
            }
            else {
                Q = H;
            }


            //R=Q.multiply(A);
        }
        assert Q != null;

        //Q=Q.transpose();
        //Q.outputToTxt("q");
       // R.outputToTxt("r");
        //QRMatrix matrix=new QRMatrix(m);
        //matrix.q=Q;
        //matrix.r=R;

        return Q.multiply(A).multiply(Q.transpose());
        //return matrix;
    }*/

    public Matrix eigenvaluesMatrixQR(double precision) throws MatrixSizeError {

       // double ak_nn;
       // assert m!=n;
       // double sum=0.;
       // double sumLastRow=0.;


       // double eps=1.e-8;

        if (m!=n){
            throw new MatrixSizeError("Матрица не квадратная");
        }
        Matrix eigenvaluesMatrix = new Matrix(m,1);
        Matrix A = this;
        int m=A.getCountRows();
        int n=A.getCountColumns();
        Matrix p;
        Matrix H;
        double sum;
        //Matrix Q=null;
        //Matrix R=A;
        int iter=0;

        //применяем метод отражений Хаусхолдера и зануляем нижние элементы матрицы лежащие ниже на 2 диагонали
        for (int k=0;k<m-2;k++){

            p=A.getPartOfMatrix(k+1,m-1,k,k);
            if (Math.abs(p.getElement(p.m-1)) <= eps)//элемент уже нулевой с определенной точностью
                continue;
            double v=p.getElement(0);//v*e_1
            if (Math.abs(v) <= eps)//элемент уже нулевой с определенной точностью
                continue;
            double norm=p.sphericalNorm();
            if (v<0)//sign(v)
                norm=-norm;
            p.setElement(0, v + norm);
            double ptp=p.transpose().multiply(p).getElement(0,0);
            H = p.multiply(p.transpose()).multiply(2./ptp);

            //H=Matrix.getEyeMatrix(m-(k+1)).subtraction(H);
            H=Matrix.getEyeMatrix(H.m).subtraction(H);
            H=Matrix.combineMatrix(getEyeMatrix(k+1),H);

            A=H.multiply(A).multiply(H.transpose());


            //R=Q.multiply(A);
        }
        do{
            //QRMatrix qrMatrix=new QRMatrix(m);
            //Matrix p;
            //Matrix H;
            Matrix Q=null;
            Matrix R=A;


            Matrix G;
            //int i_ind;
            //int j_ind;


            for (int k=0;k<m-1;k++){
                //p = R.getPartOfMatrix(0, m-1, k, k);
                //for (int j=k;j<p.m;j++) {
                    /*double v = p.getElement(j);
                    if (Math.abs(v) <= eps)//элемент уже нулевой с определенной точностью
                        continue;
                    //i=m-1 j=j
                    if (p.m-1>j) {
                        i_ind = j;
                        j_ind = p.m-1;
                    }
                    else{
                        i_ind = p.m-1;
                        j_ind = j;
                    }*/
                    int i = k;
                    int j = k+1;

                    double a_i=R.matrix[i][k];
                    double a_j=R.matrix[j][k];
                    double norm = Math.sqrt((a_i*a_i)+(a_j*a_j));
                    double cos=a_i/norm;
                    double sin=a_j/norm;

                    G=Matrix.getEyeMatrix(m);
                    G.setElement(i,i,cos);
                    G.setElement(j,j,cos);
                    G.setElement(i,j,sin);
                    G.setElement(j,i,-sin);


                    if (Q != null) {
                        Q = G.multiply(Q);
                    } else {
                        Q = G;
                    }
                    //R = Q.multiply(R);
                     R = G.multiply(R);

                   /* if (Q!=null) {
                        Q = H.multiply(Q);
                    }
                    else {
                        Q = H;
                    }*/


                    //R=Q.multiply(A);
                }
            //}
            assert Q != null;
            Q=Q.transpose();
            A=R.multiply(Q);
            //return qrMatrix;
            sum = 0.;
            //for (int i=0; i<m-1;i++)
            for (int i=1; i<m;i++)
                for (int j=0;j<i;j++)
                    sum += Math.abs(A.matrix[i][j]);

            /*for (int i=0;i<Ak.m;i++)
                Ak.matrix[i][i]+=ak_nn;
*/
        iter++;
        }
        while (sum>precision);

        //assert Q != null;

        //Q=Q.transpose();
        //Q.outputToTxt("q");
        // R.outputToTxt("r");
        //QRMatrix matrix=new QRMatrix(m);
        //matrix.q=Q;
        //matrix.r=R;


        /*for (int k=0;k<MAX_COUNT_OF_ITERATIONS;k++){
            ak_nn=Ak.matrix[Ak.m-1][Ak.n-1];
            for (int i=0;i<Ak.m;i++)
                Ak.matrix[i][i]-=ak_nn;
            QR = Ak.decompositionQR();

            Ak = QR.getR().multiply(QR.getQ());

            sumLastRow=0.;//процедура исчерпывания
            for (int i=0;i<Ak.m;i++)
                sumLastRow+=Math.abs(Ak.matrix[Ak.m-1][i]);
            if (sumLastRow<precision){
                eigenvaluesMatrix.setElement(Ak.m-1,ak_nn);
                Ak=Ak.getPartOfMatrix(0,Ak.m-2,0,Ak.n-2);
            }

            sum = 0.;
            for (int i=0; i<Ak.m-1;i++)
                for (int j=0;j<i;j++)
                    sum += Math.abs(Ak.matrix[i][j]);

            for (int i=0;i<Ak.m;i++)
                Ak.matrix[i][i]+=ak_nn;

            if (sum < precision)
                break;
        }*/

        for (int i=0;i<m;i++){
            eigenvaluesMatrix.setElement(i,A.getElement(i,i));
        }
        System.out.println("Число итераций"+iter);
        return eigenvaluesMatrix;

    }

    public LUMatrix decompositionLU() throws MatrixSizeError{

        return decompositionLU(this);
    }


    public static Matrix solveSystem(Matrix b,Matrix matrix,LUMatrix LUMatrix, double eps) throws MatrixSizeError {
        int m=matrix.getM();
        Matrix y=new Matrix(m,1);
        y.setValuesOfMatrix(0.0);
        Matrix x=generateVector(m);
        double u=0.0;
        double sum=0.0;

        y.matrix[0][0]=b.matrix[0][0];
        for (int i=1;i<m;i++) {
            sum=0.;
            for (int k = 0; k < i; k++) {
                sum += LUMatrix.getElemL(i, k) * y.matrix[k][0];
                counterSumSub++;
                counterMultDivSum++;
            }
            y.matrix[i][0] = b.matrix[i][0] - sum;
            counterSumSub++;
        }
        for (int i=m-1;i>=0;i--){
            sum=0.;
            for(int k=i+1;k<m;k++){
                sum+=LUMatrix.getElemU(i,k)*x.matrix[k][0];
                counterSumSub++;
                counterMultDivSum++;
            }
            u=LUMatrix.getElemU(i,i);
            if (Math.abs(u)>=eps) {
                x.matrix[i][0] = (y.matrix[i][0] - sum) / u;
                counterSumSub++;
                counterMultDivSum++;
            }
            else{
                System.out.println("Элемент u"+i+" "+i+"меньше eps");
                break;
            }
        }

        System.out.println("Количество операций при решении методом LU: ");
        System.out.println("'+' и '-': "+counterSumSub);
        System.out.println("'*' и '/': "+counterMultDivSum);
        int counterSum=counterMultDivSum+counterSumSub;
        System.out.println("Итого: "+counterSum);

        counterSumSub=0;
        counterMultDivSum=0;

        return x;


    }
    public Matrix solveSystem(Matrix b,LUMatrix LUMatrix) throws MatrixSizeError {
        return solveSystem(b,this,LUMatrix,this.eps);
    }



    public Matrix eigenvaluesMatrixQR_old(double precision) throws MatrixSizeError {
        int MAX_COUNT_OF_ITERATIONS=100000;
        QRMatrix QR;
        Matrix eigenvaluesMatrix = new Matrix(m,1);
        Matrix Ak = this;
        double ak_nn;
        assert m!=n;
        double sum=0.;
        double sumLastRow=0.;

        for (int k=0;k<MAX_COUNT_OF_ITERATIONS;k++){
            ak_nn=Ak.matrix[Ak.m-1][Ak.n-1];
            for (int i=0;i<Ak.m;i++)
                Ak.matrix[i][i]-=ak_nn;
            QR = Ak.decompositionQR();

            Ak = QR.getR().multiply(QR.getQ());

            sumLastRow=0.;//процедура исчерпывания
            for (int i=0;i<Ak.m;i++)
                sumLastRow+=Math.abs(Ak.matrix[Ak.m-1][i]);
            if (sumLastRow<precision){
                eigenvaluesMatrix.setElement(Ak.m-1,ak_nn);
                Ak=Ak.getPartOfMatrix(0,Ak.m-2,0,Ak.n-2);
            }

            sum = 0.;
            //for (int i=0; i<Ak.m-1;i++)
            for (int i=1; i<Ak.m;i++)
                for (int j=0;j<i;j++)
                    sum += Math.abs(Ak.matrix[i][j]);

            for (int i=0;i<Ak.m;i++)
                Ak.matrix[i][i]+=ak_nn;

            if (sum < precision)
                break;
        }

        for (int i=0;i<Ak.m;i++){
            eigenvaluesMatrix.setElement(i,Ak.getElement(i,i));
        }
        return eigenvaluesMatrix;

    }


    public void outputToScreen() {
        for (int i = 0; i < m; i++) {
            System.out.print("||| ");
            for (int j = 0; j < n; j++) {
                System.out.printf("%10.10f ", matrix[i][j]);
            }
            System.out.println(" |||\n");
        }
        System.out.println();
    }



    public void outputToTxt(String path){
        System.out.println("Writing to file...");
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
        System.out.println("Complete\n");
    }


}
