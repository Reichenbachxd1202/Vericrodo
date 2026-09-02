public class InitAuth {

    static Pairing pairing;
    static Field<Element> G0, G1, GT, Zp;

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


        Element g1 = G0.newRandomElement().getImmutable();
        Element g3 = G0.newRandomElement().getImmutable();
        Element g2 = G1.newRandomElement().getImmutable();


        for (int sig = 3; sig <= 30; sig += 3) {

            double totalTimeMs = 0.0;
            long totalApkBytes = 0;
            long totalAskBytes = 0;

            int repeatCount = 100;

            for (int repeat = 0; repeat < repeatCount; repeat++) {

                List<Element> alphaList = new ArrayList<>(sig);
                List<Element> xList = new ArrayList<>(sig);
                List<Element> AList = new ArrayList<>(sig);
                List<Element> BList = new ArrayList<>(sig);
                List<Element> CList = new ArrayList<>(sig);

                long startTime = System.nanoTime();


                for (int i = 0; i < sig; i++) {
                    alphaList.add(Zp.newRandomElement().getImmutable());
                }


                for (int i = 0; i < sig; i++) {
                    xList.add(Zp.newRandomElement().getImmutable());
                }


                Element pairingBase = pairing.pairing(g1, g2).getImmutable();


                for (int i = 0; i < sig; i++) {
                    Element A = pairingBase
                            .powZn(alphaList.get(i))
                            .getImmutable();
                    AList.add(A);
                }


                for (int i = 0; i < sig; i++) {
                    Element B = g2
                            .powZn(xList.get(i))
                            .getImmutable();
                    BList.add(B);
                }


                for (int i = 0; i < sig; i++) {
                    Element C = g1
                            .powZn(alphaList.get(i))
                            .getImmutable();
                    CList.add(C);
                }

                long endTime = System.nanoTime();
                totalTimeMs += (endTime - startTime) / 1_000_000.0;


                long apkBytes = 0;
                for (Element A : AList) {
                    apkBytes += A.toBytes().length;
                }
                for (Element B : BList) {
                    apkBytes += B.toBytes().length;
                }


                long askBytes = 0;
                for (Element C : CList) {
                    askBytes += C.toBytes().length;
                }
                for (Element x : xList) {
                    askBytes += x.toBytes().length;
                }

                totalApkBytes += apkBytes;
                totalAskBytes += askBytes;
            }

            double averageTimeMs = totalTimeMs / repeatCount;
            double averageApkBytes = (double) totalApkBytes / repeatCount;
            double averageAskBytes = (double) totalAskBytes / repeatCount;

            System.out.printf(
                    "InitAuth computation time: %.3f ms%n",
                    averageTimeMs
            );
            System.out.printf(
                    "Apk storage: %.3f Bytes%n",
                    averageApkBytes
            );
            System.out.printf(
                    "Ask storage: %.3f Bytes%n",
                    averageAskBytes
            );
            System.out.printf("sig: %d%n", sig);
            System.out.println("--------------------------------------------------");
        }


        if (g3 == null) {
            throw new IllegalStateException("Generator g3 was not initialized.");
        }
    }
}
