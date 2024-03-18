#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

// Duffing-Gleichung: y'' + δ*y' + α*y + β*y^3 = γ*cos(ω*t)

// Funktion, die das Differential der Duffing-Gleichung berechnet
void duffingEquation(double t, double y[], double dydt[]) {
    double alpha = 1.0; // Konstante α
    double beta = -1.0; // Konstante β
    double delta = 0.2; // Dämpfungskonstante δ
    double gamma = 0.3; // Amplitudenkonstante γ
    double omega = 1.0; // Frequenzkonstante ω

    dydt[0] = y[1]; // Erste Ableitung von y ist die zweite Komponente des y-Vektors
    dydt[1] = gamma * cos(omega * t) - delta * y[1] - alpha * y[0] - beta * pow(y[0], 3); // Zweite Ableitung von y
}

// Runge-Kutta-Methode vierter Ordnung zur Integration der Differentialgleichungen
void rungeKutta(double t0, double y0[], double t, double h) {
    double k1[2], k2[2], k3[2], k4[2];
    double yTemp[2];
    double dydt[2];

    ofstream output("duffing_output.txt"); // Datei zur Ausgabe öffnen

    double tNow = t0;
    double yNow[2];
    yNow[0] = y0[0];
    yNow[1] = y0[1];

    while (tNow < t) {
        // Berechne k1
        duffingEquation(tNow, yNow, dydt);
        k1[0] = h * dydt[0];
        k1[1] = h * dydt[1];

        // Berechne k2
        yTemp[0] = yNow[0] + 0.5 * k1[0];
        yTemp[1] = yNow[1] + 0.5 * k1[1];
        duffingEquation(tNow + 0.5 * h, yTemp, dydt);
        k2[0] = h * dydt[0];
        k2[1] = h * dydt[1];

        // Berechne k3
        yTemp[0] = yNow[0] + 0.5 * k2[0];
        yTemp[1] = yNow[1] + 0.5 * k2[1];
        duffingEquation(tNow + 0.5 * h, yTemp, dydt);
        k3[0] = h * dydt[0];
        k3[1] = h * dydt[1];

        // Berechne k4
        yTemp[0] = yNow[0] + k3[0];
        yTemp[1] = yNow[1] + k3[1];
        duffingEquation(tNow + h, yTemp, dydt);
        k4[0] = h * dydt[0];
        k4[1] = h * dydt[1];

        // Berechne y für den nächsten Schritt
        yNow[0] += (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6.0;
        yNow[1] += (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6.0;

        // Aktualisiere t
        tNow += h;

        // Ausgabe in Datei
        output << tNow << "\t" << yNow[0] << "\t" << yNow[1] << endl;
    }

    output.close(); // Datei schließen
}

int main() {
    double t0 = 0.0; // Anfangszeit
    double y0[2] = {0.1, 0.0}; // Anfangsbedingungen: y(0) = 0.1, y'(0) = 0
    double tFinal = 100.0; // Endzeit
    double h = 0.01; // Schrittweite

    rungeKutta(t0, y0, tFinal, h); // Aufruf der Runge-Kutta-Methode

    cout << "Simulation abgeschlossen. Ergebnisse wurden in 'duffing_output.txt' gespeichert." << endl;

    return 0;
}