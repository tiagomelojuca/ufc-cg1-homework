// ------------------------------------------------------------------------------------------------
// COMPUTACAO GRAFICA I - TAREFA 01
// ------------------------------------------------------------------------------------------------
// Esfera com o centro no eixo z negativo
// ------------------------------------------------------------------------------------------------
// Escreva um código, baseado na discussão das duas primeiras aulas,
// que pinte uma esfera no canvas. Os requisitos são os seguintes:
//   1) Defina uma janela através da qual o pintor verá a esfera 
//      >> A largura da janela, em metros, será armazenada na variável wJanela
//      >> A altura. da janela, em metros, será armazenada na variável hJanela
//      >> O centro da janela está sobre o eixo z do sistema de coordenadas
//         na posição (0, 0, - dJanela) em metros. Assim, todos os pontos no
//         plano da Janela terão coordenada z = - dJanela.
//   2) O olho do pintor está na origem do sistema de coordenadas (0,0,0)
//   3) O raio da esfera deve ser armazenado na variável rEsfera
//   4) O centro da esfera deve estar sobre o eixo z com coordenada z < - (dJanela + rEsfera)
//   5) A cor da esfera deve ser esfColor = 255, 0, 0
//   6) A cor de background deve ser cinza bgColor = 100, 100, 100
//   7) Defina o número de colunas nCol e o número de linhas nLin da matriz de cores da imagem.
//      >> nCol representa tanto o número de colunas na tela de mosquito que está presa 
//         a Janela, quanto o número de colunas do quadriculado que o pintor marcou a lápis
//         no Canvas de pintura.
//      >> nLin, de maneira análoga, representa tanto o número de linhas na tela de mosquito
//         que está presa na Janela, quanto o número de linhas do quadriculado que o pintor
//         marcou a lápis no Canvas de pintura
//      Note que
//        >> O loop aninhado mais externo é o loop de linhas que vai de 0 até nLin - 1.
//           (for int l = 0; l < nLin; l++)
//        >> O loop aninhado mais interno é o loop de colunas que vai de 0 até nCol -1
//           (for int c = 0; c < nCol; c++)
// As dimensões dos retângulos da tela de mosquito são:
//   Dx = wJanela/nCol
//   Dy = hJanela/nLin
// As coordenadas do ponto do centro de um retângulo da tela de mosquito correspondente
// ao retângulo (riscado a lápis) na posição (l, c) é dado por
//   x = - wJanela/2. + Dx/2  + c*Dx
//   y =.  hJanela/2.  -  Dy/2  -  l*Dy
// Cada raio (semi reta) parte do olho do pintor E = (0,0, 0) e  passa pelo ponto (x, y, -dJanela).
// Caso o raio enviado não tenha interseção com a esfera, armazene a cor de background na posição
// (l, c) da matriz de cores do Canvas.
// Caso o raio tangencie ou tenha interseção plena com a esfera, armazene a cor da esfera
// (esfColor) na posição (l, c) da matriz de cores do Canvas.
// ------------------------------------------------------------------------------------------------

#include <cstdint>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

// ------------------------------------------------------------------------------------------------

class TVetor4D
{
public:
    TVetor4D() = default;
    TVetor4D(double x, double y, double z, bool pt = true) : _x(x), _y(y), _z(z), _p(pt) {}

    double X() const { return _x; }
    double Y() const { return _y; }
    double Z() const { return _z; }

    void X(double x) { _x = x; }
    void Y(double y) { _y = y; }
    void Z(double z) { _z = z; }

    TVetor4D& operator+=(const TVetor4D& outro)
    {
        _x += outro._x;
        _y += outro._y;
        _z += outro._z;
        _p += outro._p;
        return *this;
    }

    TVetor4D operator+(const TVetor4D& outro) const
    {
        TVetor4D ret = *this;
        ret += outro;
        return ret;
    }

    TVetor4D& operator*=(double k)
    {
        _x *= k;
        _y *= k;
        _z *= k;
        _p *= k;
        return *this;
    }

    TVetor4D operator*(double k) const
    {
        TVetor4D ret = *this;
        ret *= k;
        return ret;
    }

    TVetor4D& operator-=(const TVetor4D& outro)
    {
        TVetor4D copiaOutro = outro;
        copiaOutro *= -1.0;
        *this += copiaOutro;
        return *this;
    }

    TVetor4D operator-(const TVetor4D& outro) const
    {
        TVetor4D ret = *this;
        ret -= outro;
        return ret;
    }

    TVetor4D& operator/=(double k)
    {
        _x /= k;
        _y /= k;
        _z /= k;
        _p /= k;
        return *this;
    }

    TVetor4D operator/(double k) const
    {
        TVetor4D ret = *this;
        ret /= k;
        return ret;
    }

protected:
    double _x = 0.0;
    double _y = 0.0;
    double _z = 0.0;
    int8_t _p =   1;
};

// ------------------------------------------------------------------------------------------------

class TPonto3D : public TVetor4D
{
public:
    TPonto3D() = default;
    TPonto3D(double x, double y, double z) : TVetor4D(x, y, z) {}
    TPonto3D(const TVetor4D& outro) : TVetor4D(outro.X(), outro.Y(), outro.Z()) {}
};

// ------------------------------------------------------------------------------------------------

class TVetor3D : public TVetor4D
{
public:
    TVetor3D() : TVetor4D(0.0, 0.0, 0.0, false) {}
    TVetor3D(double x, double y, double z) : TVetor4D(x, y, z, false) {}
    TVetor3D(const TVetor4D& outro) : TVetor4D(outro.X(), outro.Y(), outro.Z(), false) {}

    double Dot(const TVetor3D& outro) const
    {
        return _x * outro._x +
               _y * outro._y +
               _z * outro._z;
    }

    double Norma() const
    {
        return sqrt(Dot(*this));
    }

    TVetor3D Normalizado() const
    {
        TVetor3D ret = *this;
        ret /= Norma();
        return ret;
    }
};

// ------------------------------------------------------------------------------------------------

class TRaio3D
{
public:
    TRaio3D() = default;
    TRaio3D(const TPonto3D& origem, const TVetor3D& direcao) : _origem(origem), _direcao(direcao) {}

    const TPonto3D& Origem() const
    {
        return _origem;
    }
    const TVetor3D& Direcao() const
    {
        return _direcao;
    }
    TPonto3D Ponto(double t) const
    {
        return _origem + _direcao * t;
    }

private:
    TPonto3D _origem;
    TVetor3D _direcao;
};

// ------------------------------------------------------------------------------------------------

class TCor
{
public:
    TCor() = default;
    TCor(uint8_t r, uint8_t g, uint8_t b) : _r(r), _g(g), _b(b) {}

    bool operator==(const TCor& outra) const
    {
        return _r == outra._r && _g == outra._g && _b == outra._b;
    }

    uint8_t R() const { return _r; };
    uint8_t G() const { return _g; };
    uint8_t B() const { return _b; };

    void R(uint8_t r) { _r = r; };
    void G(uint8_t g) { _g = g; };
    void B(uint8_t b) { _b = b; };

private:
    uint8_t _r = 0;
    uint8_t _g = 0;
    uint8_t _b = 0;
};

// ------------------------------------------------------------------------------------------------

class IEntidade3D
{
public:
    IEntidade3D() = default;
    virtual ~IEntidade3D() = default;

    virtual IEntidade3D* Copia() const = 0;
};

// ------------------------------------------------------------------------------------------------

class TEsfera : public IEntidade3D
{
public:
    TEsfera() = delete;
    TEsfera(const TPonto3D& centro, double raio) : _centro(centro), _raio(raio) {}

    IEntidade3D* Copia() const override
    {
        return new TEsfera(*this);
    }

    const TPonto3D& Centro() const { return _centro; }
    double Raio() const { return _raio; }
    const TCor& Cor() const { return _cor; }

    void Cor(const TCor& cor) { _cor = cor; }

private:
    TPonto3D _centro;
    double _raio = 0.0;

    TCor _cor;
};

// ------------------------------------------------------------------------------------------------

class TJanela
{
public:
    TJanela() = default;
    TJanela(
        const TPonto3D& centro,
        double wJanela,
        double hJanela,
        uint16_t wCanvas,
        uint16_t hCanvas
    ) :
        _centro(centro),
        _wJanela(wJanela),
        _hJanela(hJanela),
        _wCanvas(wCanvas),
        _hCanvas(hCanvas)
    {
        _dx = _wJanela / _wCanvas;
        _dy = _hJanela / _hCanvas;
    }

    const TPonto3D& Centro() const { return _centro; }
    double Largura() const { return _wJanela; }
    double Altura() const { return _hJanela; }
    uint16_t LarguraCanvas() const { return _wCanvas; }
    uint16_t AlturaCanvas() const { return _hCanvas; }

    double X(uint16_t coluna) const
    {
        return -0.5 * _wJanela + 0.5 * _dx + coluna * _dx;
    }
    double Y(uint16_t linha) const
    {
        return 0.5 * _hJanela - 0.5 * _dy - linha * _dy;
    }

private:
    TPonto3D _centro;
    double _wJanela = 0.0;
    double _hJanela = 0.0;
    uint16_t _wCanvas = 0;
    uint16_t _hCanvas = 0;

    double _dx = 0.0;
    double _dy = 0.0;
};

// ------------------------------------------------------------------------------------------------

class TArquivoPPM
{
public:
    TArquivoPPM() = delete;
    TArquivoPPM(
        const std::string& caminho,
        uint16_t largura,
        uint16_t altura
    )
        : _caminho(caminho),
          _largura(largura),
          _altura(altura)
    {
        _stream.open(_caminho);

        if (_stream.is_open())
        {
            _stream << "P3"
                    << "\n"
                    << std::to_string(_largura)
                    << " "
                    << std::to_string(altura)
                    << "\n"
                    << "255"
                    << "\n";
        }
    }
    ~TArquivoPPM()
    {
        _stream.close();
    }

    bool Aberto() const
    {
        return _stream.is_open();
    }
    const char* Caminho() const
    {
        return _caminho.c_str();
    }

    bool Anexa(const TCor& cor)
    {
        const bool aberto = Aberto();

        if (aberto)
        {
            _stream << std::to_string(cor.R()) << " "
                    << std::to_string(cor.G()) << " "
                    << std::to_string(cor.B()) << "\n";
        }
        
        return aberto;
    }

private:
    std::string _caminho;
    uint16_t _largura;
    uint16_t _altura;

    std::ofstream _stream;
};

// ------------------------------------------------------------------------------------------------

namespace FuncoesGeometricas
{
    TVetor3D Versor(const TPonto3D& pInicio, const TPonto3D& pFim)
    {
        const TVetor3D v = pFim - pInicio;
        const TVetor3D d = v.Normalizado();

        return d;
    }

    bool Intercepta(const TRaio3D& raio, const TEsfera& esfera)
    {
        const TVetor3D w = raio.Origem() - esfera.Centro();
        const TVetor3D d = raio.Direcao();

        const double a = d.Dot(d);
        const double b = 2.0 * (w.Dot(d));
        const double c = w.Dot(w) - esfera.Raio() * esfera.Raio();
        const double delta = b * b - 4.0 * a * c;

        return delta >= 0.0;
    }
}

// ------------------------------------------------------------------------------------------------

class TCena3D
{
public:
    TCena3D() = default;
    TCena3D(const TPonto3D& origem, const TJanela& janela) : _janela(janela), _p0(origem) {}

    const TJanela& Janela() const { return _janela; }
    void Janela(const TJanela& janela) { _janela = janela; }

    const TPonto3D& Origem() const { return _p0; }
    void Origem(const TPonto3D& p0) { _p0 = p0; }

    const TCor& BgColor() const { return _bgColor; }
    void BgColor(const TCor& cor) { _bgColor = cor; }

    void Insere(const IEntidade3D& entidade)
    {
        _entidades.push_back(std::unique_ptr<IEntidade3D>(entidade.Copia()));
    }

    void Renderizar(TArquivoPPM& arq) const
    {
        const uint16_t nLinhas = _janela.AlturaCanvas();
        const uint16_t nColunas = _janela.LarguraCanvas();

        const double z = _janela.Centro().Z();
        for (int l = 0; l < nLinhas; l++)
        {
            const double y = _janela.Y(l);

            for (int c = 0; c < nColunas; c++)
            {
                const double x = _janela.X(c);

                arq.Anexa(Cor({ x, y, z }));
            }
        }
    }

private:
    TCor Cor(const TPonto3D& p) const
    {
        // essa funcao sera substituida no futuro
        // por enquanto, to deixando renderizar pela ordem de insercao dos
        // objetos na lista, e nao pela posicao relativa ao observador,
        // como deveria ser... isso esta errado, eh claro, mas deixei assim
        // enquanto nao pergunto pro prof a melhor maneira de fazer a interacao
        // entre multiplos objetos na cena
        // (talvez usar a cor do obj com menor dist da intersecao?)

        TCor pixel = _bgColor;

        for (const std::unique_ptr<IEntidade3D>& entidade : _entidades)
        {
            if (auto esfera = dynamic_cast<const TEsfera*>(entidade.get()))
            {
                if (pixel == _bgColor)
                {
                    pixel = Cor(p, *esfera);
                }
            }
        }

        return pixel;
    }

    TCor Cor(const TPonto3D& p, const TEsfera& esfera) const
    {
        const TVetor3D d = FuncoesGeometricas::Versor(_p0, p);
        const TRaio3D raio { _p0, d };
        const bool intercepta = FuncoesGeometricas::Intercepta(raio, esfera);

        return intercepta ? esfera.Cor() : _bgColor;
    }

    TJanela _janela;
    TPonto3D _p0; // olho do pintor (origem)

    TCor _bgColor;
    std::vector<std::unique_ptr<IEntidade3D>> _entidades;
};

// ------------------------------------------------------------------------------------------------

int main()
{
    const TPonto3D p0 { 0.0, 0.0, 0.0 };

    const double wJanela = 10.0;
    const double hJanela = 10.0;
    const double dJanela = 12.0;
    const uint16_t wCanvas = 500u;
    const uint16_t hCanvas = 500u;
    const TJanela janela { { 0.0, 0.0, -dJanela }, wJanela, hJanela, wCanvas, hCanvas };

    TArquivoPPM arq("teste.ppm", janela.LarguraCanvas(), janela.AlturaCanvas());

    const bool erro = !arq.Aberto();
    if (!erro)
    {
        TEsfera esfera { { 0.0, 0.0, -20.0 }, 4.0 };
        esfera.Cor({ 255u, 0u, 0u });

        // TEsfera esfera2 { { -10.0, -10.0, -55.0 }, 4.0 };
        // esfera2.Cor({ 0u, 255u, 0u });

        TCena3D cena { p0, janela };
        cena.BgColor({ 100u, 100u, 100u });
        cena.Insere(esfera);
        // cena.Insere(esfera2);

        cena.Renderizar(arq);
    }

    return erro;
}

// ------------------------------------------------------------------------------------------------
