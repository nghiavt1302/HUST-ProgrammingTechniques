#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#define INF 999.5
using namespace std;

// Tạo kiểu dữ liệu hàm đa thức
typedef struct PolynomialFunction {
    double coef[100]; // Mảng các hệ số của đơn thức bậc i
    int degree; // Bậc của đa thức
} function;

void newFileOutput(){
    ofstream newFile;
    newFile.open("output.txt", ios::trunc);

    time_t timeNow = time(0);
    char* time = ctime(&timeNow);
    newFile << time << endl;

    newFile.close();
}

function inputFunction(){
    ofstream writeToFile;
    writeToFile.open("output.txt", ios::app);

    int choice;
    function f;
    cout << "------ Input your function ------" << endl;
    cout << "1. Input function from keyboard." << endl;
    cout << "2. Input function from file." << endl;
    cout << "What do you want?" << endl;
    cin >> choice;
    switch (choice) {
        case 1:{
            cout << "Input the function's degree:";
            cin >> f.degree;
            cout << "Input the coef: " << endl;
            for (int i = 0; i < f.degree + 1; i++) {
                cout << "Coeff at x^" << i << ":";
                cin >> f.coef[i];
            }
            break;
        }
        case 2:{
            string path;
            cout << "Input the function's degree:";
            cin >> f.degree;
            ifstream inputFile;
            inputFile.open("input.txt", ios::in);
            int i = 0;
            while (!inputFile.eof()) {
                inputFile >> f.coef[i];
                i++;
            }
            inputFile.close();
            break;
        }
    }
    cout << "Your function is: ";
    writeToFile << "Your function is: ";
    for (int i = f.degree; i >=0 ; i--) {
        if (f.coef[i] < 0){
            cout << "- ";
            writeToFile << "- ";
        }
        cout << fabs(f.coef[i]) << "x^" << i << " ";
        writeToFile << fabs(f.coef[i]) << "x^" << i << " ";
        if (i==0){
            break;
        }
        if (f.coef[i-1]>=0) {
            cout << "+ ";
            writeToFile << "+ ";
        }
    }
    cout << "= 0." << endl;
    writeToFile << "= 0." << endl;
    writeToFile.close();

    return f;
}

// Tính đạo hàm
PolynomialFunction getDerivative(function f){
    function df;
    df.degree = f.degree - 1;
    for (int i = 0; i <= df.degree; i++){
        df.coef[i] = f.coef[i+1] * (i+1);
    }
    return df;
}

// Tính giá trị hàm số tại x
double funcValue(function f, double x){
    double value = 0;
    for(int i = 0; i <= f.degree; i++){
        value += f.coef[i] * pow(x, i);
    }
    return value;
}

// Tạo kiểu dữ liệu: khoảng phân ly nghiệm
typedef struct RootsRange{
    double a, b;
} rr;

// Tìm miền chứa tất cả nghiệm (định lý Horner)
rr allRootsRange(function f, double a, double b){
    rr rangeAll;
    double max_coef = fabs(f.coef[0]);
    for (int i = 1; i <= f.degree; i++) {
        if (max_coef < fabs(f.coef[i])){
            max_coef = fabs(f.coef[i]);
        }
    }

    rangeAll.a = -1 - (max_coef / fabs(f.coef[f.degree]));
    if (a > rangeAll.a){
        rangeAll.a = a;
    }

    rangeAll.b = 1 + (max_coef / fabs(f.coef[f.degree]));
    if (b < rangeAll.b){
        rangeAll.b = b;
    }

    return rangeAll;
}

// Chia miền chứa tất cả nghiệm thành các khoảng phân ly nghiệm, lưu vào mảng các khoảng phân ly nghiệm
int divideRootsRange(function f, double a, double b, rr ranges[]){
    ofstream writeToFile;
    writeToFile.open("output.txt", ios::app);

    rr rangeAll = allRootsRange(f, a, b);
    double x1 = rangeAll.a, x2 = rangeAll.b;
    function df = getDerivative(f);
    int rangeIndex = 0, count = 0;
    double temp = x1 + 0.001;
    double fx1 = funcValue(f, x1);
    double dfx1 = funcValue(df, x1), dfTemp = 0, fTemp = 0;

    while (temp <= x2){
        while (temp <= x2){
            dfTemp = funcValue(df, temp);
            fTemp = funcValue(f, temp);
            if ((dfx1*dfTemp) < 0 && (fx1*fTemp) <= 0){
                break;
            }
            temp += 0.001;
        }
        ranges[rangeIndex].a = x1;
        ranges[rangeIndex].b = temp - 0.001;
        rangeIndex++;
        count++;
        x1 = temp;
        temp += 0.001;
        dfx1 = funcValue(df, x1);
        fx1 = funcValue(f,x1);
    }

    for (int i = 0; i < rangeIndex; ++i) {
        cout << setprecision(3) << "The root range " << i + 1 << " is: [" << ranges[i].a << ", " << ranges[i].b << "]" << endl;
        writeToFile << setprecision(3) << "The root range " << i + 1 << " is: [" << ranges[i].a << ", " << ranges[i].b << "]" << endl;
    }

    cout << endl;
    writeToFile << endl;

    return count;
}

// Thu hẹp KPL nghiệm bằng phương pháp chia đôi
rr bisectionMethod(function f, int rangeIndex, rr range[], double eps){
    ofstream writeToFile;
    writeToFile.open("output.txt", ios::app);

    rr rootRangeNarrowed;
    function df = getDerivative(f), ddf = getDerivative(df);
    double a = range[rangeIndex].a;
    double b = range[rangeIndex].b;
    double x, fx, fa;
    double ddfa = funcValue(ddf, a);
    double ddfb = funcValue(ddf, b);

    while ((ddfa*ddfb < 0) || (fabs(b - a) > eps)){
        x = (a+b)/2;
        fx = funcValue(f, x);
        fa = funcValue(f, a);
        if (fa * fx > 0){
            a = x;
        } else{
            b = x;
        }
        ddfa = funcValue(ddf,a);
        ddfb = funcValue(ddf,b);
    }
    rootRangeNarrowed.a = a;
    rootRangeNarrowed.b = b;

    cout << setprecision(3) << "Bisection method: range [" << range[rangeIndex].a << ", " << range[rangeIndex].b << "]" << " is narrowed - [" << a << ", " << b << "]" << endl;
    writeToFile << setprecision(3) << "Bisection method: range [" << range[rangeIndex].a << ", " << range[rangeIndex].b << "]" << " is narrowed - [" << a << ", " << b << "]" << endl;
    
    return rootRangeNarrowed;
}

// Tính giá trị m1 = min |df(x)|
double calculateM1(function f, double a, double b){
    function df = getDerivative(f);
    double m1;
    double dfa = funcValue(df, a), dfb = funcValue(df, b);

    if (fabs(dfa) > fabs(dfb)){
        m1 = fabs(dfb);
    } else{
        m1 = fabs(dfa);
    }

    return m1;
}

// Tính giá trị M2 = max |ddf(x)|
double calculateM2(function f, double a, double b){
    double M2;
    function df = getDerivative(f), ddf = getDerivative(df);
    double ddfa = funcValue(ddf, a), ddfb = funcValue(ddf, b);

    if (fabs(ddfa) > fabs(ddfb)){
        M2 = fabs(ddfa);
    } else{
        M2 = fabs(ddfb);
    }

    return M2;
}

// Tính sai số mục tiêu
double calculateTargetError(function f, double a, double b, double x){
    double m1 = calculateM1(f, a, b);
    double targetError = fabs(funcValue(f, x))/m1;
    return targetError;
}

// Tính sai số theo hai xấp xỉ liên tiếp
double calculateTwoConsecutiveStepsError(function f, double a, double b, double x, double s){
    double m1 = calculateM1(f, a, b), M2 = calculateM2(f, a, b);
    double twoStepsError = M2 * pow(x-s,2) / (2 * m1);
    return twoStepsError;
}


void printTableHeader(int num){
    ofstream writeToFile;
    writeToFile.open("output.txt", ios::app);
    // In header bảng
    for (int i = 0; i < 15 + num + 10;i++) {
        cout << "-";
        writeToFile << "-";
    }
    cout << endl;
    writeToFile << endl;
    cout << left << setw(15) << "|  Iteration" << "|  " << left << setw(num + 6) << "x" << "|" << endl;
    writeToFile << left << setw(15) << "|  Iteration" << "|  " << left << setw(num + 6) << "x" << "|" << endl;
    for (int i = 0; i < 15 + num + 10;i++) {
        cout << "-";
        writeToFile << "-";
    }
    cout << endl;
    writeToFile << endl;
}

void printTableFooter(int num){
    ofstream writeToFile;
    writeToFile.open("output.txt", ios::app);
    // In footer bảng
    for (int i = 0; i < 15 + num + 10;i++) {
        cout << "-";
        writeToFile << "-";
    }
    cout << endl;
    writeToFile << endl;
}
// Tìm nghiệm với số lần lặp cho trước, đánh giá sai số theo 2 cách
void findRootWithConstantIteration(function f, rr rootRange, int numOfIterations, int numOfDecimals){
    ofstream writeToFile;
    writeToFile.open("output.txt", ios::app);
    function df = getDerivative(f), ddf = getDerivative(df);
    double a = rootRange.a, b = rootRange.b;
    double targetError, twoStepsError;
    int k = 0;
    double x, s;
    if (funcValue(f, a) * funcValue(ddf, a) > 0){
        x = a;
    } else{
        x = b;
    }
    printTableHeader(numOfDecimals);
    do {
        s = x;
        x = x - (funcValue(f, x) / funcValue(df, x));
        targetError = calculateTargetError(f, a, b, x);
        twoStepsError = calculateTwoConsecutiveStepsError(f, a, b, x, s);
        cout << setprecision(numOfDecimals + 1) <<  "|  " << left << setw(12)  << k + 1 << "|  " << left << setw(numOfDecimals + 6) << x << "|"<< endl;
        writeToFile << setprecision(numOfDecimals + 1) <<  "|  " << left << setw(12)  << k + 1 << "|  " << left << setw(numOfDecimals + 6) << x << "|"<< endl;
        k++;
    } while (k < numOfIterations);
    printTableFooter(numOfDecimals);
    cout << setprecision(numOfDecimals + 1) << "After " << k << " iteration, x = " << x << endl;
    cout << setprecision(4) << "Target error = " << targetError << "." << endl;
    cout << setprecision(4) << "Two consecutive steps error = " << twoStepsError << "." << endl;
    writeToFile << setprecision(numOfDecimals + 1) << "After " << k << " iteration, x = " << x << endl;
    writeToFile << setprecision(4) << "Target error = " << targetError << "." << endl;
    writeToFile << setprecision(4) << "Two consecutive steps error = " << twoStepsError << "." << endl;
}

// Tìm nghiệm theo sai số e cho trước theo công thức sai số mục tiêu
void findRootWithGivenError1(function f, rr rootRange, double error, int numOfDecimals){
    ofstream writeToFile;
    writeToFile.open("output.txt", ios::app);
    function df = getDerivative(f), ddf = getDerivative(df);
    double a = rootRange.a, b = rootRange.b;
    int k = 0;
    double s, x;
    double m1 = calculateM1(f, a, b);
    if (funcValue(f, a) * funcValue(ddf, a) > 0){
        x = a;
    } else{
        x = b;
    }
    printTableHeader(numOfDecimals);
    do {
        s = x;
        x = x - (funcValue(f, x) / funcValue(df, x));
//        cout << setprecision(numOfDecimals + 1) << "Iteration " << k + 1 << ", x = " << x << endl;
        cout << setprecision(numOfDecimals + 1) <<  "|  " << left << setw(12)  << k + 1 << "|  " << left << setw(numOfDecimals + 6) << x << "|"<< endl;
        writeToFile << setprecision(numOfDecimals + 1) <<  "|  " << left << setw(12)  << k + 1 << "|  " << left << setw(numOfDecimals + 6) << x << "|"<< endl;
        k++;
    } while ((m1 * error) <= fabs(funcValue(f, x)));
    printTableFooter(numOfDecimals);
    cout << setprecision(numOfDecimals + 1) << "With target error = " << error << ", after " << k << " iteration, x = " << x << endl;
    writeToFile << setprecision(numOfDecimals + 1) << "With target error = " << error << ", after " << k << " iteration, x = " << x << endl;
}

// Tìm  theo nghiệm sai số e cho trước theo công thức sai số hai xấp xỉ liên tiếp
void findRootWithGivenError2(function f, rr rootRange, double error, int numOfDecimals){
    ofstream writeToFile;
    writeToFile.open("output.txt", ios::app);
    function df = getDerivative(f), ddf = getDerivative(df);
    double a = rootRange.a, b = rootRange.b;
    int k = 0;
    double s, x;
    double m1 = calculateM1(f, a, b), M2 = calculateM2(f, a, b);
    if (funcValue(f, a) * funcValue(ddf, a) > 0){
        x = a;
    } else{
        x = b;
    }
    printTableHeader(numOfDecimals);
    double err = sqrt(fabs(2.0 * m1 * error / M2));
    do {
        s = x;
        x = x - (funcValue(f, x) / funcValue(df, x));
//        cout << setprecision(numOfDecimals + 1) << "Iteration " << k + 1 << ", x = " << x << endl;
        cout << setprecision(numOfDecimals + 1) <<  "|  " << left << setw(12)  << k + 1 << "|  " << left << setw(numOfDecimals + 6) << x << "|"<< endl;
        writeToFile << setprecision(numOfDecimals + 1) <<  "|  " << left << setw(12)  << k + 1 << "|  " << left << setw(numOfDecimals + 6) << x << "|"<< endl;
        k++;
    } while (fabs(x - s) >= err);
    printTableFooter(numOfDecimals);
    cout << setprecision(numOfDecimals + 1) << "With two consecutive steps error = " << error << ", after " << k << " iteration, x = " << x << endl;
    writeToFile << setprecision(numOfDecimals + 1) << "With two consecutive steps error = " << error << ", after " << k << " iteration, x = " << x << endl;
}

// Tìm nghiệm với sai số e cho trước thoả mãn |x(n) - x(n-1)| <= e
void findRootWithGivenErrorBetweenToSteps(function f, rr rootRange, double error, int numOfDecimals){
    ofstream writeToFile;
    writeToFile.open("output.txt", ios::app);
    function df = getDerivative(f), ddf = getDerivative(df);
    double a = rootRange.a, b = rootRange.b;
    int k = 0;
    double s, x;
    if (funcValue(f, a) * funcValue(ddf, a) > 0){
        x = a;
    } else{
        x = b;
    }
    printTableHeader(numOfDecimals);
    do {
        s = x;
        x = x - (funcValue(f, x) / funcValue(df, x));
//        cout << setprecision(numOfDecimals + 1) << "Iteration " << k + 1 << ", x = " << x << endl;
        cout << setprecision(numOfDecimals + 1) <<  "|  " << left << setw(12)  << k + 1 << "|  " << left << setw(numOfDecimals + 6) << x << "|"<< endl;
        writeToFile << setprecision(numOfDecimals + 1) <<  "|  " << left << setw(12)  << k + 1 << "|  " << left << setw(numOfDecimals + 6) << x << "|"<< endl;
        k++;
    } while (fabs(x - s) > error);
    printTableFooter(numOfDecimals);
    cout << setprecision(numOfDecimals + 1) << "With two consecutive steps error = " << error << ", after " << k << " iteration, x = " << x << endl;
    writeToFile << setprecision(numOfDecimals + 1) << "With two consecutive steps error = " << error << ", after " << k << " iteration, x = " << x << endl;
}

void readFileOutput(){
    ifstream readFile;
    readFile.open("output.txt", ios::in);
    string  string1;
    for(istreambuf_iterator<char,char_traits<char> > it(readFile.rdbuf());
        it != istreambuf_iterator<char,char_traits<char> >(); it++)   {

        string1 += *it; //append at end of string
    }
    cout << string1.data() << endl;
    readFile.close();
}

int main() {
    newFileOutput();
    cout << "----- Welcome -----" << endl;
    ofstream writeToFile;
    writeToFile.open("output.txt", ios::app);
    function f = inputFunction();
    int choice;
    int endLoop;
    rr range[f.degree];
    int numOfSolutions = divideRootsRange(f, -INF, INF, range);
    rr rangeNarrow[numOfSolutions];
    for (int i = 0; i < numOfSolutions; i++) {
        rangeNarrow[i] = bisectionMethod(f, i, range, 0.5);
    }
    cout << endl;
    writeToFile << endl;
    do {
        cout << "1. Find roots with given number of iterations." << endl;
        cout << "2. Find roots with given error by target error." << endl;
        cout << "3. Find roots with given error by two consecutive steps error." << endl;
        cout << "4. Find roots with given error where |x(n) - x(n-1)| < error." << endl;
        cin >> choice;
        switch (choice) {
            case 1:{
                writeToFile << "1. Find roots with given number of iterations." << endl;
                int numOfIteration, numDecimal;
                cout << "How many iterations:";
                cin >> numOfIteration;
                writeToFile << "Iterations: " << numOfIteration << "." << endl;
                cout << "How many decimals:";
                cin >> numDecimal;
                writeToFile << "Number of decimals: " << numDecimal << "." << endl;
                for (int i = 0; i < numOfSolutions; i++) {
                    cout << "----- Solution " << i + 1 << " -----" << endl;
                    writeToFile << "----- Solution " << i + 1 << " -----" << endl;
                    findRootWithConstantIteration(f, rangeNarrow[i], numOfIteration, numDecimal);
                }
                cout << endl;
                writeToFile << endl;
                cout << "Press 1 to redo. Press any to stop" << endl;
                cin >> endLoop;
                break;
            }
            case 2:{
                writeToFile << "2. Find roots with given error by target error." << endl;
                int numDecimal;
                double err;
                cout << "Error:";
                cin >> err;
                writeToFile << "Error: " << err << endl;
                cout << "How many decimals:";
                cin >> numDecimal;
                writeToFile << "Number of decimals: " << numDecimal << "." << endl;
                for (int i = 0; i < numOfSolutions; i++) {
                    cout << "----- Solution " << i + 1 << " -----" << endl;
                    writeToFile << "----- Solution " << i + 1 << " -----" << endl;
                    findRootWithGivenError1(f, rangeNarrow[i], err, numDecimal);
                }
                cout << endl;
                writeToFile << endl;
                cout << "Press 1 to redo. Press any to stop" << endl;
                cin >> endLoop;
                break;
            }
            case 3:{
                writeToFile << "3. Find roots with given error by two consecutive steps error." << endl;
                int numDecimal;
                double err;
                cout << "Error:";
                cin >> err;
                writeToFile << "Error: " << err << endl;
                cout << "How many decimals:";
                cin >> numDecimal;
                writeToFile << "Number of decimals: " << numDecimal << "." << endl;
                for (int i = 0; i < numOfSolutions; i++) {
                    cout << "----- Solution " << i + 1 << " -----" << endl;
                    writeToFile << "----- Solution " << i + 1 << " -----" << endl;
                    findRootWithGivenError2(f, rangeNarrow[i], err, numDecimal);
                }
                cout << endl;
                writeToFile << endl;
                cout << "Press 1 to redo. Press any to stop" << endl;
                cin >> endLoop;
                break;
            }
            case 4:{
                writeToFile << "4. Find roots with given error where |x(n) - x(n-1)| < error." << endl;
                int numDecimal;
                double err;
                cout << "Error:";
                cin >> err;
                writeToFile << "Error: " << err << endl;
                cout << "How many decimals:";
                cin >> numDecimal;
                writeToFile << "Number of decimals: " << numDecimal << "." << endl;
                for (int i = 0; i < numOfSolutions; i++) {
                    cout << "----- Solution " << i + 1 << " -----" << endl;
                    writeToFile << "----- Solution " << i + 1 << " -----" << endl;
                    findRootWithGivenErrorBetweenToSteps(f, rangeNarrow[i], err, numDecimal);
                }
                cout << endl;
                writeToFile << endl;
                cout << "Press 1 to redo. Press any to stop" << endl;
                cin >> endLoop;
                break;
            }
        }
    } while (endLoop == 1);
    writeToFile.close();
    cout << "See your result in output.txt." << endl;
    int seeResult;
    cout << "Press 1 to read results from output.txt file.";
    cin >> seeResult;
    if (seeResult == 1){
        system("cls");
        readFileOutput();
    }
}
