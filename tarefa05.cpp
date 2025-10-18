// ------------------------------------------------------------------------------------------------
// COMPUTACAO GRAFICA I - TAREFA 05
// ------------------------------------------------------------------------------------------------
// Esta cena contem os seguintes objetos:
// >> cinco planos
// >> um cilindro
// >> um cone
// >> um cubo 
// >> uma esfera
// Alem disso, contem uma luz pontual e uma luz ambiente.
// Objetivo: Gerar a cena colorida como resultado das interacoes entre as propriedades dos
//           materiais e as fontes luminosas.
// Plano 1: Chao
// >> Ponto P_pi = (0, -150cm, 0)
// >> Vetor unitario normal ao plano: n = (0, 1., 0.)
// >> Kd = Ke = Ka = Textura de madeira

// Plano 2: Parede lateral direita 
// >> Ponto P_pi = (200cm, -150cm, 0)
// >> Vetor unitario normal ao plano: n = (-1., 0., 0.)
// >> Kd = Ke = Ka = (0.686, 0.933, 0.933)

// Plano 3: Parede frontal 
// >> Ponto P_pi = (200cm, -150cm, -400cm)
// >> Vetor unitario normal ao plano: n = (0., 0., 1.)
// >> Kd = Ke = Ka = (0.686, 0.933, 0.933)

// Plano 4: Parede lateral esquerda
// >> Ponto P_pi = (-200cm, -150cm, 0cm)
// >> Vetor unitario normal ao plano: n = (1., 0., 0.)
// >> Kd = Ke = Ka = (0.686, 0.933, 0.933)

// Plano 5: Teto
// >> Ponto P_pi = (0cm, 150cm, 0cm)
// >> Vetor unitario normal ao plano: n = (0., -1., 0.)
// >> Kd = Ke = Ka = (0.933, 0.933, 0.933)

// Cilindro:
// >> Centro da base: C_b = (0cm, -150cm, -200cm)
// >> Raio do cilindro: R = 5cm
// >> Altura do cilindro: H = 90cm
// >> Vetor unitario do eixo do cilindro: n = (0., 1., 0.)
// >> Kd = Ke = Ka = (0.824, 0.706, 0.549)

// Cone:
// >> Centro da base: C_b = (0cm, -60cm, -200cm)
// >> Raio da base do cone: R = 90cm
// >> Altura do cone: H = 150cm
// >> Vetor unitario do eixo do cone: n = (0., 1., 0.)
// >> Kd = Ke = Ka = (0., 1., 0.498)

// Cubo:
// >> Aresta do cubo: a = 40cm
// >> Centro da base do cubo: C_c = (0cm, -150cm, -165cm)
// >> As arestas da base sao paralelas aos eixos x e z do sistema de coordenadas
// >> Kd = Ke = Ka = (1., 0.078., 0.576)


// Esfera:
// >> Centro da esfera: C = (0, 95cm, -200)
// >> Raio da esfera: R = 5cm
// >> Kd = Ke = Ka = (0.854, 0.647, 0.125)

// Fonte de luz pontual:
// >> Intensidade da fonte:  I_F = (0.7, 0.7, 0.7)
// >> Posicao da fonte:  P_F = (-100cm, 140cm, -20cm)

// Fonte de luz ambiente>> Intensidade da fonte:  I_F = (0.3, 0.3, 0.3)


// Novidades:  Textura do chao e escala e translacao do cubo.
// ------------------------------------------------------------------------------------------------
