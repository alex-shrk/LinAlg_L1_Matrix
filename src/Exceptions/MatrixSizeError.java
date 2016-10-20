package Exceptions;

public class MatrixSizeError extends Exception {
    public MatrixSizeError() {
        super("Несовместимые размеры матриц");
    }
    public MatrixSizeError(String str){
        super(str);
    }
}
