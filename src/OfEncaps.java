public class OfEncaps {

    static Pairing pairing;
    static Field<Element> G0, G1, GT, Zp;


    public static Element hashFromBytesToG0(byte[] bytes) {
        return G0.newElement()
                .setFromHash(bytes, 0, bytes.length)
                .getImmutable();
    }


    public static int hashFromUToTheta(Element element, int[] Theta) {
        if (Theta == null || Theta.length == 0) {
            throw new IllegalArgumentException("Theta must not be empty.");
        }

        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            byte[] hash = digest.digest(element.toBytes());
            int index = new BigInteger(1, hash)
                    .mod(BigInteger.valueOf(Theta.length))
                    .intValue();
            return Theta[index];
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException("SHA-256 is not available.", e);
        }
    }

    public static void main(String[] args) {


        List<String> curveList = Arrays.asList(
                "g149.properties",
                "a.properties"
        );

        for (String curveFile : curveList) {
            System.out.printf("Curve: %s%n", curveFile);
            runForCurve(curveFile);
            System.out.println("==================================================");
        }
    }

    private static void runForCurve(String curveFile) {

        pairing = PairingFactory.getPairing(curveFile);


        G0 = pairing.getG1();
        G1 = pairing.getG2();
        GT = pairing.getGT();
        Zp = pairing.getZr();

        Random random = new Random();


        Element g1 = G0.newRandomElement().getImmutable();
        Element g3 = G0.newRandomElement().getImmutable();
        Element g2 = G1.newRandomElement().getImmutable();


        Element M = GT.newRandomElement().getImmutable();

        int repeat = 30;


        for (int thetaSize = 3; thetaSize <= 30; thetaSize += 3) {


            int[] Theta = new int[thetaSize];
            for (int sig = 1; sig <= thetaSize; sig++) {
                Theta[sig - 1] = sig;
            }

            Element[] alpha = new Element[thetaSize];
            Element[] x = new Element[thetaSize];
            Element[] A = new Element[thetaSize];
            Element[] B = new Element[thetaSize];
            Element[] C = new Element[thetaSize];


            for (int sig = 1; sig <= thetaSize; sig++) {
                alpha[sig - 1] = Zp.newRandomElement().getImmutable();
            }


            for (int sig = 1; sig <= thetaSize; sig++) {
                x[sig - 1] = Zp.newRandomElement().getImmutable();
            }


            Element authorityPairingBase = pairing
                    .pairing(g1, g2)
                    .getImmutable();


            for (int sig = 1; sig <= thetaSize; sig++) {
                A[sig - 1] = authorityPairingBase
                        .powZn(alpha[sig - 1])
                        .getImmutable();
            }


            for (int sig = 1; sig <= thetaSize; sig++) {
                B[sig - 1] = g2
                        .powZn(x[sig - 1])
                        .getImmutable();
            }


            for (int sig = 1; sig <= thetaSize; sig++) {
                C[sig - 1] = g1
                        .powZn(alpha[sig - 1])
                        .getImmutable();
            }


            for (int N = 10; N <= 100; N += 10) {

                double totalTimeMs = 0.0;
                long totalICBytes = 0;

                for (int rep = 0; rep < repeat; rep++) {


                    Element v = Zp.newRandomElement().getImmutable();
                    Element[] w = new Element[N];
                    Element[] xi = new Element[N];
                    Element[] t = new Element[N];

                    for (int i = 1; i <= N; i++) {
                        w[i - 1] = Zp.newRandomElement().getImmutable();
                    }
                    for (int i = 1; i <= N; i++) {
                        xi[i - 1] = Zp.newRandomElement().getImmutable();
                    }
                    for (int i = 1; i <= N; i++) {
                        t[i - 1] = Zp.newRandomElement().getImmutable();
                    }

                    Element[][] tldC_i_sig_2 = new Element[N][thetaSize];
                    Element[] tldC_i_3 = new Element[N];
                    Element[][] tldC_i_sig_4 = new Element[N][thetaSize];
                    Element[] tldC_i_5 = new Element[N];

                    long startTime = System.nanoTime();


                    Element pairingBase = pairing
                            .pairing(g1, g2)
                            .getImmutable();


                    Element tldC1 = M
                            .duplicate()
                            .mul(pairingBase.powZn(v))
                            .getImmutable();


                    for (int i = 1; i <= N; i++) {
                        Element e_w_i = pairingBase
                                .powZn(w[i - 1])
                                .getImmutable();

                        for (int sig = 1; sig <= thetaSize; sig++) {
                            Element A_sig_t_i = A[sig - 1]
                                    .powZn(t[i - 1])
                                    .getImmutable();
                            tldC_i_sig_2[i - 1][sig - 1] = e_w_i
                                    .duplicate()
                                    .mul(A_sig_t_i)
                                    .getImmutable();
                        }
                    }


                    for (int i = 1; i <= N; i++) {
                        tldC_i_3[i - 1] = g2
                                .powZn(t[i - 1])
                                .getImmutable();
                    }


                    for (int i = 1; i <= N; i++) {
                        for (int sig = 1; sig <= thetaSize; sig++) {
                            tldC_i_sig_4[i - 1][sig - 1] = B[sig - 1]
                                    .powZn(t[i - 1])
                                    .getImmutable();
                        }
                    }


                    for (int i = 1; i <= N; i++) {
                        tldC_i_5[i - 1] = g2
                                .powZn(xi[i - 1])
                                .getImmutable();
                    }

                    long endTime = System.nanoTime();
                    totalTimeMs += (endTime - startTime) / 1_000_000.0;



                    long icBytes = tldC1.toBytes().length;

                    for (int i = 1; i <= N; i++) {
                        for (int sig = 1; sig <= thetaSize; sig++) {
                            icBytes += tldC_i_sig_2[i - 1][sig - 1]
                                    .toBytes().length;
                        }

                        icBytes += tldC_i_3[i - 1].toBytes().length;

                        for (int sig = 1; sig <= thetaSize; sig++) {
                            icBytes += tldC_i_sig_4[i - 1][sig - 1]
                                    .toBytes().length;
                        }

                        icBytes += tldC_i_5[i - 1].toBytes().length;
                    }

                    totalICBytes += icBytes;
                }

                double averageTimeMs = totalTimeMs / repeat;
                double averageICKBytes = (double) totalICBytes
                        / repeat
                        / 1024.0;

                System.out.printf("Theta=%d, N=%d%n", thetaSize, N);
                System.out.printf(
                        "OfEncaps computation time: %.3f ms%n",
                        averageTimeMs
                );
                System.out.printf(
                        "IC storage: %.3f KBytes%n",
                        averageICKBytes
                );
                System.out.println("--------------------------------------------------");
            }


            if (Theta.length != thetaSize
                    || A.length != thetaSize
                    || B.length != thetaSize
                    || C.length != thetaSize) {
                throw new IllegalStateException("Authority initialization failed.");
            }
        }


        if (g3 == null || random == null) {
            throw new IllegalStateException("Global initialization failed.");
        }
    }
}
