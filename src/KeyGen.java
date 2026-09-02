public class KeyGen {

    static Pairing pairing;
    static Field<Element> G0, G1, GT, Zp;


    public static Element hashToG0(byte[] data) {
        return G0.newElement()
                .setFromHash(data, 0, data.length)
                .getImmutable();
    }


    public static int hashToTheta(Element u, List<Integer> Theta) {
        if (Theta == null || Theta.isEmpty()) {
            throw new IllegalArgumentException("Theta must not be empty.");
        }

        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            byte[] hash = digest.digest(u.toBytes());
            int index = new BigInteger(1, hash)
                    .mod(BigInteger.valueOf(Theta.size()))
                    .intValue();
            return Theta.get(index);
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

        SecureRandom random = new SecureRandom();


        Element g1 = G0.newRandomElement().getImmutable();
        Element g3 = G0.newRandomElement().getImmutable();
        Element g2 = G1.newRandomElement().getImmutable();


        byte[] gid = new byte[32];
        random.nextBytes(gid);

        int repeat = 30;


        for (int sig = 3; sig <= 30; sig += 3) {

            Map<Integer, Element> alphaMap = new HashMap<>(sig);
            Map<Integer, Element> xMap = new HashMap<>(sig);
            Map<Integer, Element> AMap = new HashMap<>(sig);
            Map<Integer, Element> BMap = new HashMap<>(sig);
            Map<Integer, Element> CMap = new HashMap<>(sig);


            for (int j = 1; j <= sig; j++) {
                alphaMap.put(j, Zp.newRandomElement().getImmutable());
            }


            for (int j = 1; j <= sig; j++) {
                xMap.put(j, Zp.newRandomElement().getImmutable());
            }


            Element pairingBase = pairing.pairing(g1, g2).getImmutable();


            for (int j = 1; j <= sig; j++) {
                Element A_j = pairingBase
                        .powZn(alphaMap.get(j))
                        .getImmutable();
                Element B_j = g2
                        .powZn(xMap.get(j))
                        .getImmutable();
                Element C_j = g1
                        .powZn(alphaMap.get(j))
                        .getImmutable();

                AMap.put(j, A_j);
                BMap.put(j, B_j);
                CMap.put(j, C_j);
            }


            for (int uni = 10; uni <= 100; uni += 10) {

                double totalTimeMs = 0.0;
                long totalAtkBytes = 0;
                long totalCekBytes = 0;

                for (int rep = 0; rep < repeat; rep++) {


                    List<byte[]> Delta = new ArrayList<>(uni);
                    for (int i = 1; i <= uni; i++) {
                        byte[] delta_i = new byte[32];
                        random.nextBytes(delta_i);
                        Delta.add(delta_i);
                    }

                    Map<Integer, Element> rMap = new HashMap<>(sig);
                    Map<Integer, Element> atk_j_1_Map = new HashMap<>(sig);
                    Map<Integer, Element> atk_j_2_Map = new HashMap<>(sig);
                    Map<Integer, Map<Integer, Element>> atk_i_j_3_Map =
                            new HashMap<>(uni);
                    Map<Integer, Element> tauMap = new HashMap<>(sig);
                    Map<Integer, Element> cek_j_1_Map = new HashMap<>(sig);
                    Map<Integer, Element> cek_j_2_Map = new HashMap<>(sig);

                    long startTime = System.nanoTime();


                    for (int j = 1; j <= sig; j++) {
                        rMap.put(j, Zp.newRandomElement().getImmutable());
                    }


                    for (int j = 1; j <= sig; j++) {
                        Element H_gid = hashToG0(gid);
                        Element H_gid_x_j = H_gid
                                .powZn(xMap.get(j))
                                .getImmutable();
                        Element atk_j_1 = CMap.get(j)
                                .duplicate()
                                .mul(H_gid_x_j)
                                .getImmutable();
                        atk_j_1_Map.put(j, atk_j_1);
                    }


                    for (int j = 1; j <= sig; j++) {
                        Element atk_j_2 = BMap.get(j)
                                .powZn(rMap.get(j))
                                .getImmutable();
                        atk_j_2_Map.put(j, atk_j_2);
                    }


                    for (int i = 1; i <= uni; i++) {
                        Map<Integer, Element> atk_j_3_Map = new HashMap<>(sig);
                        for (int j = 1; j <= sig; j++) {
                            Element H_delta_i = hashToG0(Delta.get(i - 1));
                            Element atk_i_j_3 = H_delta_i
                                    .powZn(rMap.get(j))
                                    .getImmutable();
                            atk_j_3_Map.put(j, atk_i_j_3);
                        }
                        atk_i_j_3_Map.put(i, atk_j_3_Map);
                    }


                    for (int j = 1; j <= sig; j++) {
                        tauMap.put(j, Zp.newRandomElement().getImmutable());
                    }


                    for (int j = 1; j <= sig; j++) {
                        Element cek_j_1 = atk_j_1_Map.get(j)
                                .powZn(tauMap.get(j))
                                .getImmutable();
                        cek_j_1_Map.put(j, cek_j_1);
                    }


                    for (int j = 1; j <= sig; j++) {
                        Element cek_j_2 = g2
                                .powZn(tauMap.get(j))
                                .getImmutable();
                        cek_j_2_Map.put(j, cek_j_2);
                    }

                    long endTime = System.nanoTime();
                    totalTimeMs += (endTime - startTime) / 1_000_000.0;


                    long atkBytes = 0;
                    for (Element element : atk_j_1_Map.values()) {
                        atkBytes += element.toBytes().length;
                    }
                    for (Element element : atk_j_2_Map.values()) {
                        atkBytes += element.toBytes().length;
                    }
                    for (Map<Integer, Element> atk_j_3_Map
                            : atk_i_j_3_Map.values()) {
                        for (Element element : atk_j_3_Map.values()) {
                            atkBytes += element.toBytes().length;
                        }
                    }


                    long cekBytes = 0;
                    for (Element element : cek_j_1_Map.values()) {
                        cekBytes += element.toBytes().length;
                    }
                    for (Element element : cek_j_2_Map.values()) {
                        cekBytes += element.toBytes().length;
                    }

                    totalAtkBytes += atkBytes;
                    totalCekBytes += cekBytes;
                }

                double averageTimeMs = totalTimeMs / repeat;
                double averageAtkBytes = (double) totalAtkBytes / repeat;
                double averageCekBytes = (double) totalCekBytes / repeat;

                System.out.printf("sig=%d, uni=%d%n", sig, uni);
                System.out.printf(
                        "KeyGen computation time: %.3f ms%n",
                        averageTimeMs
                );
                System.out.printf(
                        "Atk storage: %.3f Bytes%n",
                        averageAtkBytes
                );
                System.out.printf(
                        "Cek storage: %.3f Bytes%n",
                        averageCekBytes
                );
                System.out.println("--------------------------------------------------");
            }


            if (AMap.size() != sig) {
                throw new IllegalStateException("Authority initialization failed.");
            }
        }


        if (g3 == null) {
            throw new IllegalStateException("Generator g3 was not initialized.");
        }
    }
}
