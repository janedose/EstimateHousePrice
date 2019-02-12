import java.util.*;
import java.io.*;
/**
 * This program predicts the price of a house based on historical data in an input file.
 */
public class EstimateHousePrice {
    public static final String INPUTFILE = "historicaldata.txt";
    /**
     * The main method to calculate parameters based on historical data, get user input, 
     * and predict the price of a house.
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public static void main(String[] args) 
            throws FileNotFoundException, IOException {
        Scanner console = new Scanner (System.in);
        Scanner input = new Scanner(new File(INPUTFILE));
        int m = rowCount(input);
        int[][]table=dataTable(input, m);
        double[]normalizedSize = normalize(m, table, 0);
        double[]normalizedBeds = normalize(m, table, 1);
        double[][] X= X(m,normalizedSize,normalizedBeds);
        double[][] Y= Y(m,table);
        double[][]theta = theta(m,X,Y);
        predict(console, theta, table, m);
        System.out.println("The calculated parameters are: " + Arrays.deepToString(theta));
    }
    /**
     * This method counts the number of rows in the input file.
     * @param input The data from input file
     * @return the number of rows in the input file
     * @throws IOException 
     */
    public static int rowCount(Scanner input) throws IOException {
        LineNumberReader lnr = new LineNumberReader(new FileReader(INPUTFILE));
        int numRows=0;
        while (lnr.readLine()!=null) {
            numRows++;
        }
        lnr.close();
        return numRows;
    }
    /**
     * This method organizes data from the input file into a table.
     * @param input The data from input file
     * @param numRows the number of rows in the input file
     * @return input data in a table format
     * @throws IOException 
     */
    public static int[][] dataTable(Scanner input, int numRows) throws IOException {
        int[][]table = new int [numRows][3];
        int row = 0;
        while (input.hasNextLine()) {
        String line = input.nextLine();
        String[] lineData = line.split(",");
        for (int i=0; i<3; i++) {
            table[row][i]=Integer.parseInt(lineData[i]);
        }
        row++;
        }
        return table;
    } 
    /**
     * This method calculates mu for normalizing data.
     * @param numHouses the number of houses
     * @param table input data in a table format
     * @param j column index
     * @return mu to normalize data
     */
    public static double mu(int numHouses, int[][] table, int j) {
        double mu = 0;
        for (int i = 0; i <numHouses; i++) {
           mu += table[i][j];
        }
        mu=mu/numHouses;
        return mu;
    }
    /**
     * This method calculates sigma for normalizing data.
     * @param numHouses the number of houses
     * @param table input data in a table format
     * @param j column index
     * @param mu the calculated mu to normalize data
     * @return sigma to normalize data
     */
    public static double Sigma(int numHouses, int[][] table, int j, double mu) {
        double smallSigma = 0;
        for (int i = 0; i <numHouses; i++) {
           smallSigma += (table[i][j]-mu)*(table[i][j]-mu);
        }
        smallSigma=Math.sqrt(smallSigma/(numHouses-1));
        return smallSigma;
    }
    /**
     * This method normalizes the data from the input file.
     * @param numHouses the number of houses
     * @param table input data in a table format
     * @param j column index
     * @return the normalized data in a table format
     */
    public static double[] normalize(int numHouses, int[][] table, int j) {
        double mu = mu(numHouses, table, j);      
        double smallSigma = Sigma(numHouses, table, j, mu);
        double[]normalized=new double[numHouses];
        for (int i = 0; i <numHouses; i++) {
            normalized[i]=(table[i][j]-mu)/smallSigma;
        }
        return normalized;
    }
    /**
     * This method derives the matrix X from the normalized data.
     * @param m number of houses
     * @param normalizedSize the normalized size of houses
     * @param normalizedBeds the normalized number of beds of houses
     * @return the matrix X
     */
    public static double[][] X(int m, double[] normalizedSize, double[] normalizedBeds) {
        double[][]matrixX = new double[m][3];
        for (int i=0; i<m; i++) {
            matrixX[i][0]=1;
        }
        for (int i=0; i<m; i++) {
            matrixX[i][1]=normalizedSize[i];
        }
        for (int i=0; i<m; i++) {
            matrixX[i][2]=normalizedBeds[i];
        }
        return matrixX;
    }
    /**
     * This method derives the matrix Y from the input file for house prices.
     * @param m number of houses
     * @param table input data in a table format
     * @return the matrix Y
     */
    public static double[][] Y(int m, int[][] table) {
        double[][]matrixY = new double[m][1];
        for (int i=0; i<m; i++) {
            matrixY[i][0]=table[i][2];
        }
        return matrixY;
    }
    /**
     * This method transposes the matrix.
     * @param matrix the matrix to be transposed
     * @return the resulting matrix
     */
    public static double[][] transpose(double [][] matrix){
        double[][] TMatrix = new double[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix.length; i++)
            for (int j = 0; j < matrix[0].length; j++)
                TMatrix[j][i] = matrix[i][j];
        return TMatrix;
    }
    /**
     * This method multiplies 2 matrices.
     * @param matrix1 the first matrix
     * @param matrix2 the second matrix
     * @return the resulting matrix
     */
    public static double[][] multiply(double[][] matrix1, double[][] matrix2) { 
        if(matrix1[0].length != matrix2.length) return null;
        double[][] MMatrix = new double[matrix1.length][matrix2[0].length];
        for(int i = 0; i < matrix1.length; i++) {         
            for(int j = 0; j < matrix2[0].length; j++) {     
                for(int k = 0; k < matrix1[0].length; k++) { 
                    MMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
        return MMatrix;
    }
    /**
     * This method subtracts 2 matrices.
     * @param matrix1 the first matrix
     * @param matrix2 the second matrix
     * @return the resulting matrix
     */
    public static double[][] subtract(double[][] matrix1, double[][] matrix2) {
        if (matrix1.length != matrix2.length || matrix1[0].length != matrix2[0].length) return null;
        double[][] SMatrix = new double[matrix1.length][matrix1[0].length];
        for(int i=0;i<matrix1.length;i++) {
            for(int j=0;j<matrix1[0].length;j++) {
                SMatrix[i][j] = matrix1[i][j] - matrix2[i][j];
            }
        }
        return SMatrix;
    }
    /**
     * This method calculates the parameters for predicting house prices.
     * @param m the number of houses
     * @param X the matrix X
     * @param Y the matrix Y
     * @return the parameters
     */
    public static double[][] theta(int m, double[][] X, double[][] Y) {
        double[][]thetaOld=new double[3][1];
        double[][]theta=new double[3][1];
        double[][]thetaTemp=new double[3][1];
        double[][]XT=transpose(X);
        double[][]YT=transpose(Y);
        for (int times=0;times<400;times++) {
            thetaTemp = transpose(multiply(subtract(multiply(transpose(thetaOld),XT),YT),X));
            for (int i=0;i<3;i++) {
            theta[i][0] = thetaOld[i][0]-thetaTemp[i][0]*.01/m;
            }
            thetaOld=theta;
        }
        return theta;
    }
    /**
     * This method predicts the price of a house with given size and number of bedrooms.
     * @param console The Scanner object to retrieve user inputs
     * @param theta the parameters for predicting house prices
     * @param table input data in a table format
     * @param m the number of houses
     */
    public static void predict(Scanner console, double[][] theta, int[][]table, int m) {
        System.out.println("Enter house size in square feet:");
        while (!console.hasNextDouble()) {
            console.next();
            System.out.println("Invalid input. Enter house size in square feet:");
        }
        double size = console.nextDouble();
        System.out.println("Enter number of bedrooms:");
        while (!console.hasNextInt()) {
            console.next();
            System.out.println("Invalid input. Enter number of bedrooms in integer format:");
        }
        int beds = console.nextInt();
        double muSize = mu(m, table, 0);
        double muBeds = mu(m, table, 1);
        double sigmaSize = Sigma(m, table, 0, muSize);
        double sigmaBeds = Sigma(m, table, 1, muBeds);
        double price = theta[0][0]+theta[1][0]*((size-muSize)/sigmaSize)+theta[2][0]*((beds-muBeds)/sigmaBeds);
        price = Math.round(price*100.0)/100.0;
        System.out.println("\nThe predicted price is $" +price);
    }
} 