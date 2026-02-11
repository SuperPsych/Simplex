import java.util.*;

public class Main {
    public static void main(String[] args) {
        double[][] A = {
                {1, 1, 2},
                {2, 0, 3},
                {2, 1, 3}
        };
        double[] b = {4, 5, 7};
        double[] c = {3, 2, 4};

        Simplex simplex = new Simplex(A, b, c);
        double res = simplex.solve();
        System.out.println(res);
    }
}

class Simplex{
    Matrix matrix;
    Set<Integer> basics = new HashSet<>();
    public Simplex(double[][] m) {
        matrix = new Matrix(m);
    }
    public Simplex(double[][] A, double[] b, double[] c){
        this(constructMatrix(A,b,c));
    }
    public static double[][] constructMatrix(double[][] A, double[] b, double[] c){
        double[][] m = new double[A.length+1][A[0].length+A.length+2];
        m[0][0] = 1;
        for(int i = 0; i < A.length; i++){
            for(int j = 0; j < A[0].length; j++){
                m[i+1][j+1] = A[i][j];
            }
        }
        for(int i = 0; i < A.length; i++){
            m[i+1][i+A[0].length+1] = 1;
            m[i+1][m[0].length-1] = b[i];
        }
        for(int i = 0; i < c.length; i++){
            m[0][i+1] = -c[i];
        }
        return m;
    }
    public boolean nullIsFeasible(){
        for(int row=1; row<matrix.rows; row++){
            if(matrix.vals[row][matrix.cols-1]<0) return false;
        }
        return true;
    }
    public int[] getArtificialRows(){
        int n = 0;
        for(int row=1; row<matrix.rows; row++){
            if(matrix.vals[row][matrix.cols-1]<0) n++;
        }
        int[] rows = new int[n];
        int i=0;
        for(int row=1; row<matrix.rows; row++){
            if(matrix.vals[row][matrix.cols-1]<0) rows[i++] = row;
        }
        return rows;
    }
    public Set<Integer> findFeasibleBasis(){
        Simplex artificialSimplex = constructArtificialSimplex();
        double res = artificialSimplex.compute();
        if(res!=0) return null;
        return artificialSimplex.basics;
    }
    public Simplex constructArtificialSimplex(){
        int[] artificialRows = getArtificialRows();
        double[][] artificialMatrix = new double[matrix.rows][matrix.cols+artificialRows.length];
        artificialMatrix[0][0] = 1;
        for(int row=1; row<matrix.rows; row++){
            for(int col=0; col<matrix.cols-1; col++){
                artificialMatrix[row][col] = matrix.vals[row][col];
            }
        }
        int i=matrix.cols-1;
        for(int row : artificialRows){
            artificialMatrix[row][i++] = -1;
        }
        for(int row=1; row<matrix.rows; row++){
            artificialMatrix[row][artificialMatrix[0].length-1] = matrix.vals[row][matrix.cols-1];
        }
        for(i = matrix.cols-1; i<matrix.cols+artificialRows.length-1; i++){
            artificialMatrix[0][i] = 1;
        }
        Simplex artificialSimplex = new Simplex(artificialMatrix);
        for(int row : artificialRows){
            artificialSimplex.matrix.rowAdd(1,row,0);
            artificialSimplex.matrix.rowScale(-1,row);
        }
        i=0;
        for(int row=1; row<matrix.rows; row++){
            if(i<artificialRows.length && artificialRows[i]==row){
                artificialSimplex.basics.add(matrix.cols-1+i);
                i++;
            }
            else{
                artificialSimplex.basics.add(matrix.cols-matrix.rows-2+row);
            }
        }
        return artificialSimplex;
    }
    public void normalizeMatrix(){
        int[][] pivots = new int[basics.size()][2];
        int i=0;
        for(int basic : basics){
            pivots[i][0] = i+1;
            pivots[i][1] = basic;
            i++;
        }
        matrix.pivotify(pivots);
    }
    public double solve(){
        if(nullIsFeasible()){
            for(int i=matrix.cols-matrix.rows; i<matrix.cols-1; i++){
                basics.add(i);
            }
            return compute();
        }
        basics = findFeasibleBasis();
        if(basics == null) return Double.MIN_VALUE;
        normalizeMatrix();
        return compute();
    }
    public double compute(){
        int pivot = choosePivot();
        while(pivot != -1) {
            int[] kick = findKickPivot(pivot);
            int kickRow = kick[0];
            int kickCol = kick[1];
            if (kickRow == -1) {
                return Double.MIN_VALUE;
            }
            basics.remove(kickCol);
            basics.add(pivot);
            int[][] pivots = {{kickRow, pivot}};
            matrix.pivotify(pivots);
            pivot = choosePivot();
        }
        return matrix.vals[0][matrix.cols-1];
    }
    public int choosePivot(){
        for(int i=1; i<matrix.cols-1; i++){
            if(matrix.vals[0][i]<0){
                return i;
            }
        }
        return -1;
    }
    public int[] findKickPivot(int pivot){
        double min = Double.MAX_VALUE;
        int kickRow = -1;
        for(int row=1; row<matrix.rows; row++){
            double ratio = matrix.vals[row][matrix.cols-1]/matrix.vals[row][pivot];
            if(ratio < 0) continue;
            if(ratio < min){
                min = ratio;
                kickRow = row;
            }
        }
        int kickCol = -1;
        for(int basic : basics){
            if(matrix.vals[kickRow][basic]==1){
                kickCol = basic;
                break;
            }
        }
        return new int[]{kickRow, kickCol};
    }
}

class Matrix{
    double[][] vals;
    int cols;
    int rows;
    public Matrix(double[][] v){
        vals = v.clone();
        cols = v[0].length;
        rows = v.length;
    }
    public Matrix(int a, int b){
        cols = a;
        rows = b;
        vals = new double[rows][cols];
    }
    public void pivotify(int[][] pivots){
        for(int[] pivot : pivots){
            clear(pivot);
            normalize(pivot);
        }
    }
    public void clear(int[] pivot){
        for(int row=0; row<rows; row++){
            if(pivot[0] == row) continue;
            double scalar = -vals[row][pivot[1]]/vals[pivot[0]][pivot[1]];
            rowAdd(scalar, pivot[0], row);
        }
    }
    public void normalize(int[] pivot){
        rowScale(1/vals[pivot[0]][pivot[1]], pivot[0]);
    }
    public void rowAdd(double scalar, int src, int dest){
        for(int col=0; col<cols; col++){
            vals[dest][col] += scalar * vals[src][col];
        }
    }
    public void rowScale(double scalar, int row){
        for(int col=0; col<cols; col++){
            vals[row][col] *= scalar;
        }
    }
    public String toString(){
        StringBuilder sb = new StringBuilder();
        for(int row=0; row<rows; row++){
            for(int col=0; col<cols; col++){
                sb.append(vals[row][col]).append(" ");
            }
            sb.append("\n");
        }
        return sb.toString();
    }
}
