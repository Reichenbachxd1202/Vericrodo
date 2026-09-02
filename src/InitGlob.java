public class InitGlob {

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
            System.out.println("===========================================");
        }
    }

    private static void runForCurve(String curveFile) {

        pairing = PairingFactory.getPairing(curveFile);


        G0 = pairing.getG1();
        G1 = pairing.getG2();
        GT = pairing.getGT();
        Zp = pairing.getZr();

        Random random = new Random();
        int repeat = 100;


        for (int u = 10; u <= 100; u += 10) {


            for (int t = 3; t <= 30; t += 3) {

                double totalDuration = 0.0;

                for (int r = 0; r < repeat; r++) {


                    long startTime = System.nanoTime();


                    Element g1 = G0.newRandomElement().getImmutable();
                    Element g3 = G0.newRandomElement().getImmutable();
                    Element g2 = G1.newRandomElement().getImmutable();


                    Element pairingResult = pairing.pairing(g1, g2).getImmutable();


                    Element[] U = new Element[u];
                    for (int i = 0; i < u; i++) {
                        U[i] = Zp.newRandomElement().getImmutable();
                    }


                    int[] Theta = new int[t];
                    for (int i = 0; i < t; i++) {
                        Theta[i] = i + 1;
                    }


                    byte[] random256 = new byte[32];
                    random.nextBytes(random256);
                    Element hValue = hashFromBytesToG0(random256);


                    int[] fValues = new int[u];
                    for (int i = 0; i < u; i++) {
                        fValues[i] = hashFromUToTheta(U[i], Theta);
                    }

                    long endTime = System.nanoTime();
                    totalDuration += (endTime - startTime) / 1_000_000.0;


                    if (g3 == null || pairingResult == null || hValue == null
                            || fValues.length != u) {
                        throw new IllegalStateException("InitGlob failed.");
                    }
                }

                double averageDuration = totalDuration / repeat;
                System.out.printf("u = %d, t = %d%n", u, t);
                System.out.printf(
                        "InitGlob computation time: %.2f ms%n",
                        averageDuration
                );
                System.out.println("-------------------------------------------");
            }
        }
    }
}
