import java.util.Scanner;

public class ChlorideInitiation {

    private static final double SECONDS_PER_YEAR = 365.25 * 24 * 3600;
    private static final double TREF = 28 * 24 * 3600;
    private static final double PI = Math.PI;

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        System.out.print("Choose surface chloride model:\n" +
                         "1. Constant surface concentration\n" +
                         "2. Linearly increasing (f(t) = k·t)\n" +
                         "3. Sub-linearly increasing (f(t) = k·√t)\n" +
                         "Enter choice (1/2/3): ");
        int model = scanner.nextInt();

        System.out.print("Water-cement ratio (w/c): ");
        double wc = scanner.nextDouble();

        System.out.print("Fly ash percentage (%): ");
        double fa = scanner.nextDouble();

        System.out.print("Ground Granulated Blast furnace Slag Percentage (%): ");
        double sg = scanner.nextDouble();

        System.out.print("Critical chloride threshold Ct (% wt.): ");
        double ct = scanner.nextDouble();

        System.out.print("Cover depth (mm): ");
        double x = scanner.nextDouble() / 1000.0; 

        System.out.print("Diffusion coefficient at 28 days D28 (m²/s): ");
        double d28 = scanner.nextDouble();

        System.out.print("Maxm Surface chloride concentration: ");
        double csMax = scanner.nextDouble();

        System.out.print("Time to reach Maxm Surface chloride concentration for getting k value (years): ");
        double tMax = scanner.nextDouble();

        // Compute surface function constant k
        double k = csMax / (
            model == 2 ? tMax * SECONDS_PER_YEAR :
            model == 3 ? Math.sqrt(tMax * SECONDS_PER_YEAR) :
            1.0);

        // Time exponent m from fly ash + slag formula
        double m = 0.2 + 0.4 * ((fa / 50.0) + (sg / 70.0));

        double initiationTime;
        if (model == 1) {
            initiationTime = calculateConstantModel(d28, m, csMax, ct, x);
        } else if (model == 2) {
            initiationTime = calculateLinearModel(d28, m, k, ct, x);
        } else if (model == 3) {
            initiationTime = calculateSubLinearModel(d28, m, k, ct, x);
        } else {
            System.out.println("Invalid model selected.");
            scanner.close();
            return;
        }

        System.out.println("\n===== FINAL RESULT =====");
        System.out.printf("Time exponent (m): %.3f\n", m);
        System.out.printf("Chloride initiation time: %.2f years\n", initiationTime);
        if (initiationTime >= 2000) {
            System.out.println("Warning: Initiation time is at or beyond the search limit (2000 years). Actual initiation time may be longer.");
        }

        scanner.close();
    }

    // Model 1: Constant Cs (classic erfc)
    private static double calculateConstantModel(double d28, double m, double cs, double ct, double x) {
        double low = 0.1 * SECONDS_PER_YEAR;
        double high = 2000 * SECONDS_PER_YEAR;

        for (int i = 0; i < 300; i++) {
            double t = 0.5 * (low + high);
            double dt = d28 * Math.pow(TREF / t, m);
            double c = cs * erfc(x / Math.sqrt(4 * dt * t));

            System.out.printf("Iter %3d: t=%.2f years, c=%.6f, Ct=%.6f\n", i + 1, t / SECONDS_PER_YEAR, c, ct);

            if (Math.abs(c - ct) < 1e-6)
                return t / SECONDS_PER_YEAR;
            if (c < ct)
                low = t;
            else
                high = t;
        }
        return 0.5 * (low + high) / SECONDS_PER_YEAR;
    }

    // Model 2: Cs = k·t
    private static double calculateLinearModel(double d28, double m, double k, double ct, double x) {
        double low = 0.1 * SECONDS_PER_YEAR;
        double high = 2000 * SECONDS_PER_YEAR;

        for (int i = 0; i < 300; i++) {
            double t = 0.5 * (low + high);
            double dt = d28 * Math.pow(TREF / t, m);
            double sqrtDtT = Math.sqrt(4 * dt * t);
            double erfcTerm = erfc(x / sqrtDtT);
            double expTerm = Math.exp(-x * x / (4 * dt * t));
            double part1 = (1 + (x * x) / (2 * dt * t)) * erfcTerm;
            double part2 = (x / Math.sqrt(PI * dt * t)) * expTerm;
            double c = k * t * (part1 - part2);

            System.out.printf("Iter %3d: t=%.2f years, c=%.6f, Ct=%.6f\n", i + 1, t / SECONDS_PER_YEAR, c, ct);

            if (Math.abs(c - ct) < 1e-6)
                return t / SECONDS_PER_YEAR;
            if (c < ct)
                low = t;
            else
                high = t;
        }
        return 0.5 * (low + high) / SECONDS_PER_YEAR;
    }

    // Model 3: Cs = k·√t
    private static double calculateSubLinearModel(double d28, double m, double k, double ct, double x) {
        double low = 0.1 * SECONDS_PER_YEAR;
        double high = 2000 * SECONDS_PER_YEAR;

        for (int i = 0; i < 300; i++) {
            double t = 0.5 * (low + high);
            double dt = d28 * Math.pow(TREF / t, m);
            double sqrtDtT = Math.sqrt(4 * dt * t);
            double arg = x / sqrtDtT;
            double c = k * Math.sqrt(t) * (
                        Math.exp(-x * x / (4 * dt * t)) -
                        (x / Math.sqrt(4 * dt * t)) * erfc(arg));

            System.out.printf("Iter %3d: t=%.2f years, c=%.6f, Ct=%.6f\n", i + 1, t / SECONDS_PER_YEAR, c, ct);

            if (Math.abs(c - ct) < 1e-6)
                return t / SECONDS_PER_YEAR;
            if (c < ct)
                low = t;
            else
                high = t;
        }
        return 0.5 * (low + high) / SECONDS_PER_YEAR;
    }

    // Error function complement (erfc)
    private static double erfc(double z) {
        double t = 1.0 / (1.0 + 0.5 * Math.abs(z));
        double polynomial = t * (1.00002368 +
                t * (0.37409196 +
                t * (0.09678418 +
                t * (-0.18628806 +
                t * (0.27886807 +
                t * (-1.13520398 +
                t * (1.48851587 +
                t * (-0.82215223 +
                t * 0.17087277))))))));
        double ans = t * Math.exp(-z * z - 1.26551223 + polynomial);
        return (z >= 0) ? ans : 2.0 - ans;
    }
}
