// ------------------------------------------------------------------------------------------------
// COMPUTACAO GRAFICA I - TAREFA 02
// ------------------------------------------------------------------------------------------------
// Esfera iluminada por uma fonte de luz pontual
// ------------------------------------------------------------------------------------------------
// Modificar o metodo da Tarefa 01 para que, caso haja intersecao de um raio com a esfera,
// a cor retornada seja dada pela energia luminosa que vem do ponto de intersecao PI em direcao
// ao olho do observador. Essa energia luminosa e o resultado da interacao entre a energia
// luminosa emitida pela fonte pontual e o material da esfera no ponto de intersecao.
// 
// Ela eh composta de duas parcelas: a reflexao DIFUSA (I_d) e a reflexao Especular. (I_e), onde
// I_d =( I_F@K)* (l . n)
// I_e = (I_F@K)*(v . r)^m.
// 
// Use os seguintes atributos da fonte de luz pontual:
// I_F = (0.7, 0.7, 0.7)  // Intensidade da fonte pontual
// P_F = (0, 5, 0)   // Posição da fonte pontual situada a 5 metros acima do olho do observador.
// ------------------------------------------------------------------------------------------------

#include <cstdint>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

#include "BitmapPlusPlus.hpp"

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

    void Clamp(double min, double max)
    {
        ClampX(min, max);
        ClampY(min, max);
        ClampZ(min, max);
    }
    void ClampX(double min, double max)
    {
        _x = std::clamp(_x, min, max);
    }
    void ClampY(double min, double max)
    {
        _y = std::clamp(_y, min, max);
    }
    void ClampZ(double min, double max)
    {
        _z = std::clamp(_z, min, max);
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

    TVetor3D Arroba(const TVetor3D& outro) const
    {
        return {
            _x * outro._x,
            _y * outro._y,
            _z * outro._z
        };
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

class TMaterial
{
public:
    TMaterial() = default;

    double KdR() const
    {
        return _kdR;
    }
    void KdR(double r)
    {
        _kdR = r;
    }
    double KdG() const
    {
        return _kdG;
    }
    void KdG(double g)
    {
        _kdG = g;
    }
    double KdB() const
    {
        return _kdB;
    }
    void KdB(double b)
    {
        _kdB = b;
    }

    double KeR() const
    {
        return _keR;
    }
    void KeR(double r)
    {
        _keR = r;
    }
    double KeG() const
    {
        return _keG;
    }
    void KeG(double g)
    {
        _keG = g;
    }
    double KeB() const
    {
        return _keB;
    }
    void KeB(double b)
    {
        _keB = b;
    }

private:
    double _kdR = 0.0;
    double _kdG = 0.0;
    double _kdB = 0.0;

    double _keR = 0.0;
    double _keG = 0.0;
    double _keB = 0.0;
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

    const TMaterial& Material() const { return _material; }
    void Material(const TMaterial& material) { _material = material; }

private:
    TPonto3D _centro;
    double _raio = 0.0;

    TMaterial _material;
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

enum class EFormatoImagem { LOG, PPM, BMP };

// ------------------------------------------------------------------------------------------------

std::string Extensao(EFormatoImagem formato)
{
    if (formato == EFormatoImagem::LOG)
    {
        return "log";
    }

    if (formato == EFormatoImagem::PPM)
    {
        return "ppm";
    }

    if (formato == EFormatoImagem::BMP)
    {
        return "bmp";
    }

    return "";
}

// ------------------------------------------------------------------------------------------------

class IArquivoSaida
{
public:
    IArquivoSaida() = default;
    virtual ~IArquivoSaida() = default;

    virtual bool Aberto() const = 0;
    virtual const char* Caminho() const = 0;
    virtual bool Anexa(const TCor& cor) = 0;
    virtual void Flush() = 0;
};

// ------------------------------------------------------------------------------------------------

class TArquivoLOG : public IArquivoSaida
{
public:
    TArquivoLOG() = delete;
    TArquivoLOG(
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
            _stream << "W = " << std::to_string(_largura)
                    << "; "
                    << "H = "
                    << std::to_string(altura)
                    << ";\n";
        }
    }
    ~TArquivoLOG()
    {
        _stream.close();
    }

    bool Aberto() const override
    {
        return _stream.is_open();
    }
    const char* Caminho() const override
    {
        return _caminho.c_str();
    }

    bool Anexa(const TCor& cor) override
    {
        const bool aberto = Aberto();

        if (aberto)
        {
            _stream << "X = " << std::to_string(_coordX) << "; "
                    << "Y = " << std::to_string(_coordY) << "; "
                    << "R = " << std::to_string(cor.R()) << "; "
                    << "G = " << std::to_string(cor.G()) << "; "
                    << "B = " << std::to_string(cor.B()) << ";\n";
            
            _coordX++;
            if (_coordX == _largura)
            {
                _coordX = 0;
                _coordY++;
            }
        }
        
        return aberto;
    }

    bool Anexa(const std::string& logMsg)
    {
        const bool aberto = Aberto();

        if (aberto)
        {
            _stream << "    " << logMsg << "\n";
        }
        
        return aberto;
    }

    void Flush() override
    {
        _stream.flush();
    }

private:
    std::string _caminho;
    uint16_t _largura;
    uint16_t _altura;

    std::ofstream _stream;
    uint16_t _coordX = 0;
    uint16_t _coordY = 0;
};

// ------------------------------------------------------------------------------------------------

class TArquivoPPM : public IArquivoSaida
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

    bool Aberto() const override
    {
        return _stream.is_open();
    }
    const char* Caminho() const override
    {
        return _caminho.c_str();
    }

    bool Anexa(const TCor& cor) override
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

    void Flush() override
    {
        _stream.flush();
    }

private:
    std::string _caminho;
    uint16_t _largura;
    uint16_t _altura;

    std::ofstream _stream;
};

// ------------------------------------------------------------------------------------------------

class TArquivoBMP : public IArquivoSaida
{
public:
    TArquivoBMP() = delete;
    TArquivoBMP(
        const std::string& caminho,
        uint16_t largura,
        uint16_t altura
    )
        : _caminho(caminho),
          _largura(largura),
          _altura(altura)
    {
        _imagem = new bmp::Bitmap(_largura, _altura);
    }
    ~TArquivoBMP()
    {
        delete _imagem;
    }

    bool Aberto() const override
    {
        return _imagem->width() > 0 && _imagem->height() > 0;
    }
    const char* Caminho() const override
    {
        return _caminho.c_str();
    }

    bool Anexa(const TCor& cor) override
    {
        const bool aberto = Aberto();

        if (aberto)
        {
            _imagem->set(x, y, bmp::Pixel { cor.R(), cor.G(), cor.B() });

            x++;
            if (x == _imagem->width())
            {
                x = 0;
                y++;
            }
        }
        
        return aberto;
    }

    void Flush() override
    {
        _imagem->save(_caminho);
    }

private:
    std::string _caminho;
    uint16_t _largura;
    uint16_t _altura;

    bmp::Bitmap* _imagem = nullptr;
    std::int32_t x = 0;
    std::int32_t y = 0;
};

// ------------------------------------------------------------------------------------------------

namespace FuncoesGerais
{
    uint8_t Trunca(double n)
    {
        return static_cast<uint8_t>(floor(n));
    }

    TCor Vec2Cor(const TVetor3D& v)
    {
        return { Trunca(v.X()), Trunca(v.Y()), Trunca(v.Z()) };
    }

    std::unique_ptr<IArquivoSaida> FabricaArquivo(
        EFormatoImagem formato,
        const std::string& nomeArquivo,
        uint16_t largura,
        uint16_t altura
    )
    {
        IArquivoSaida* arq = nullptr;

        std::string _nomeArquivo = nomeArquivo;
        _nomeArquivo += ".";
        _nomeArquivo += Extensao(formato);

        if (formato == EFormatoImagem::LOG)
        {
            arq = new TArquivoLOG(_nomeArquivo, largura, altura);
        }
        else if (formato == EFormatoImagem::PPM)
        {
            arq = new TArquivoPPM(_nomeArquivo, largura, altura);
        }
        else if (formato == EFormatoImagem::BMP)
        {
            arq = new TArquivoBMP(_nomeArquivo, largura, altura);
        }

        return std::unique_ptr<IArquivoSaida>(arq);
    }
}

// ------------------------------------------------------------------------------------------------

namespace FuncoesGeometricas
{
    TVetor3D Versor(const TPonto3D& pInicio, const TPonto3D& pFim)
    {
        const TVetor3D v = pFim - pInicio;
        const TVetor3D d = v.Normalizado();

        return d;
    }

    std::vector<double> Intersecoes(const TRaio3D& raio, const TEsfera& esfera)
    {
        std::vector<double> intersecoes;

        const TVetor3D w = raio.Origem() - esfera.Centro();
        const TVetor3D d = raio.Direcao();

        const double a = d.Dot(d);
        const double b = 2.0 * (w.Dot(d));
        const double c = w.Dot(w) - esfera.Raio() * esfera.Raio();
        const double delta = b * b - 4.0 * a * c;

        if (delta < 0.0)
        {
            return intersecoes;
        }

        intersecoes.push_back((-b - sqrt(delta)) / 2.0 * a);
        intersecoes.push_back((-b + sqrt(delta)) / 2.0 * a);

        return intersecoes;
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

    void Renderizar(IArquivoSaida& arq)
    {
        if (auto arqLog = dynamic_cast<TArquivoLOG*>(&arq))
        {
            _arqLog = arqLog;
        }

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

        _arqLog = nullptr;
        arq.Flush();
    }

    void Log(const std::string& msg) const
    {
        if (_arqLog != nullptr)
        {
            _arqLog->Anexa(msg);
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
        const TMaterial& material = esfera.Material();
        const std::vector<double> intersecoes = FuncoesGeometricas::Intersecoes(raio, esfera);

        if (intersecoes.empty())
        {
            return _bgColor;
        }

        const TPonto3D pFonte = { 0.0, 5.0, 0.0 };
        const TPonto3D p1 = raio.Ponto(intersecoes[0]);

        const TVetor3D n = TVetor3D { p1 - esfera.Centro() } / esfera.Raio();
        const TVetor3D l = TVetor3D { pFonte - p1 }.Normalizado();
        const TVetor3D v = d * -1.0;
        const TVetor3D r = (n * 2.0 * n.Dot(l)) - l;

        const TPonto3D iFonte { 0.7, 0.7, 0.7 };
        const TPonto3D kd { material.KdR(), material.KdG(), material.KdB() };
        const TPonto3D ke { material.KeR(), material.KeG(), material.KeB() };

        const double m = 10.0;
        const double fd = std::max(0.0, n.Dot(l));
        const double fe = pow(std::max(v.Dot(r), 0.0), m);

        const auto id = TVetor3D(iFonte).Arroba(kd) * fd;
        const auto ie = TVetor3D(iFonte).Arroba(ke) * fe;
        const TVetor3D i = id + ie;
        TVetor3D c = i.Arroba(TVetor3D(255.0, 255.0, 255.0));
        c.Clamp(0.0, 255.0);

        return FuncoesGerais::Vec2Cor(c);
    }

    TArquivoLOG* _arqLog = nullptr;

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

    std::unique_ptr<IArquivoSaida> arq = FuncoesGerais::FabricaArquivo(
        EFormatoImagem::BMP, "teste", janela.LarguraCanvas(), janela.AlturaCanvas()
    );

    const bool erro = !arq->Aberto();
    if (!erro)
    {
        TMaterial material;
        material.KdR(1.0);
        material.KdG(0.0);
        material.KdB(0.0);
        material.KeR(1.0);
        material.KeG(1.0);
        material.KeB(1.0);

        TEsfera esfera { { 0.0, 0.0, -20.0 }, 4.0 };
        esfera.Material(material);

        // TEsfera esfera2 { { -10.0, -10.0, -55.0 }, 4.0 };
        // esfera2.Cor({ 0u, 255u, 0u });

        TCena3D cena { p0, janela };
        cena.BgColor({ 100u, 100u, 100u });
        cena.Insere(esfera);

        // cena.Insere(esfera2);

        cena.Renderizar(*arq);
    }

    return erro;
}

// ------------------------------------------------------------------------------------------------
