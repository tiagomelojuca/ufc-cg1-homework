// ------------------------------------------------------------------------------------------------
// COMPUTACAO GRAFICA I - TRABALHO FINAL
// ------------------------------------------------------------------------------------------------
// 1.     Definicao do cenario 
// 1.1.    Coerencia tematica (Obrigatorio)
// O cenario tem de ter coerencia tematica, isto e, nao pode ser um amontoado de objetos
// aleatoriamente distribuidos.
// 1.2.    Coordenadas do mundo (Obrigatorio) 
// O cenario deve ser montado de forma que todos os objetos estejam no primeiro octante, isto e,
// as coordenadas dos vertices de todos os objetos terao x, y e z positivos.
// 1.3.    Objetos
// 1.3.1. Tipos de objetos (apresentar pelo menos um objeto de cada tipo) (Obrigatorio)
// ·       Esfera
// ·       Cilindro
// ·       Cone
// ·       Malha
// 1.3.2. Materiais (pelo menos quatro materiais distintos) (Obrigatorio)
// 1.3.3. Textura (pelo menos uma textura aplicada) (Obrigatorio)
// 1.4.    Transformacoes
// 1.4.1. Translacao (Obrigatorio)
// 1.4.2. Rotacao
// ·       Em torno de um dos eixos x, y ou z (Obrigatorio)
// ·       Em torno de um eixo arbitrario (Obrigatorio um dos dois metodos)
//              o    Mudanca de sistemas de coordenadas 
//              o    Quaternios
// 1.4.3. Escala (Obrigatorio)
// 1.4.4. Cisalhamento (+ 0.5)
// 1.4.5. Espelho em relacao a um plano arbitrario (+ 0.5)
// 1.5.    Fontes luminosas
// 1.5.1. Pontual (Obrigatorio)
// 1.5.2. Spot (+1.0)
// 1.5.3. Direcional (+0.5)
// 1.5.4. Ambiente (Obrigatorio)
// 2.     Camera
// 2.1.    Permitir a especificacao de (Obrigatorio)
// 2.1.1. Posicao da camera (Eye)
// 2.1.2. Direcionamento de visada (At point)
// 2.1.3. Orientacao da camara em torno do eixo de visada (Up point)
// 2.2.    Parametros adicionais (Obrigatorio)
// 2.2.1. Distancia focal (d)
// 2.2.2. Campo de visao (definir as coordenadas de camera da janela: xmin, xmax, ymin, ymax) 
// 3.     Projecoes
// 3.1.    Perspectiva (Obrigatorio)
// 3.1.1. Alterar os parametros adicionais da camera para  
// 3.1.1.1. aumentar o campo de visao (zoom out) (Obrigatorio)
// 3.1.1.2. diminuir o campo de visao (zoom in)(Obrigatorio)
// 3.1.2. Demonstrar como posicionar a camera para obter
// 3.1.2.1. Perspectiva com um ponto de fuga (+0.5)
// 3.1.2.2. Perspectiva com dois pontos de fuga (+ 0.5)
// 3.1.2.3. Perspectiva com tres ou mais pontos de fuga (+ 0.5)
// 3.2.    Ortografica (+ 0.5)
// 3.3.    Obliqua (+0.5)
// 4.     Sombra (Obrigatorio)
// 5.     Interatividade
// 5.1.    Implementar a funcao de pick (Obrigatorio)
// 5.2.    Uso de interface grafica (Bonus de 0.5 a 1.0)
// 6.     Imagem gerada por ray casting com pelo menos 500 x 500 pixels (Obrigatorio)
// 7.     Bonus de criatividade e beleza (ate 1.0)
// ------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <memory>
#include <utility>

#ifndef UNICODE
#define UNICODE
#endif
#include <windows.h>
#include <windowsx.h>

#include "BitmapPlusPlus.hpp"

// ------------------------------------------------------------------------------------------------
static bool bp = false;
// ------------------------------------------------------------------------------------------------

template <typename T, typename S>
class TMatriz
{
public:
    class TLinha // Smart Handle
    {
    public:
        TLinha(T* ptrLinha) : _ptr(ptrLinha) {}

        T& operator[](S coluna)
        {
            return _ptr[coluna - 1];
        }

    private:
        T* _ptr = nullptr;
    };

    TMatriz(S linhas, S colunas)
    {
        _linhas = linhas;
        _colunas = colunas;
        AlocaMem(_linhas, _colunas);
    }

    TMatriz(const TMatriz& outra)
    {
        _linhas = outra._linhas;
        _colunas = outra._colunas;
        AlocaMem(_linhas, _colunas);

        for (S linha = 0; linha < _linhas; linha++)
        {
            for (S coluna = 0; coluna < _colunas; coluna++)
            {
                _matriz[linha][coluna] = outra._matriz[linha][coluna];
            }
        }
    }

    TMatriz(TMatriz&& outra)
    {
        _linhas = outra._linhas;
        outra._linhas = 0;

        _colunas = outra._colunas;
        outra._colunas = 0;

        _matriz = outra._matriz;
        outra._matriz = nullptr;
    }

    TMatriz& operator=(const TMatriz& outra)
    {
        if (&outra != this)
        {
            _linhas = outra._linhas;
            _colunas = outra._colunas;
            AlocaMem(_linhas, _colunas);

            for (S linha = 0; linha < _linhas; linha++)
            {
                for (S coluna = 0; coluna < _colunas; coluna++)
                {
                    _matriz[linha][coluna] = outra._matriz[linha][coluna];
                }
            }
        }

        return *this;
    }

    TMatriz& operator=(TMatriz&& outra)
    {
        if (&outra != this)
        {
            _linhas = outra._linhas;
            outra._linhas = 0;

            _colunas = outra._colunas;
            outra._colunas = 0;

            _matriz = outra._matriz;
            outra._matriz = nullptr;
        }

        return *this;
    }

    ~TMatriz()
    {
        for (S linha = 0; linha < _linhas; linha++)
        {
            delete[] _matriz[linha];
        }

        delete[] _matriz;
    }

    void Inicializa(const T& valorPadrao)
    {
        for (S linha = 0; linha < _linhas; linha++)
        {
            for (S coluna = 0; coluna < _colunas; coluna++)
            {
                _matriz[linha][coluna] = valorPadrao;
            }
        }
    }

    TLinha operator[](S linha)
    {
        return TLinha(_matriz[linha - 1]);
    }

    S NumeroLinhas() const
    {
        return _linhas;
    }
    S NumeroColunas() const
    {
        return _colunas;
    }

protected:
    void AlocaMem(S linhas, S colunas)
    {
        _matriz = new T*[linhas];
        for (S linha = 0; linha < linhas; linha++)
        {
            _matriz[linha] = new T[colunas];
        }
    }

    S _linhas;
    S _colunas;

    T** _matriz;
};

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

    TVetor3D Vetorial(const TVetor3D& outro) const
    {
        return TVetor3D{
            _y * outro._z - _z * outro._y,
            _z * outro._x - _x * outro._z,
            _x * outro._y - _y * outro._x
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

class IDispositivoSaida
{
public:
    IDispositivoSaida() = default;
    virtual ~IDispositivoSaida() = default;

    virtual bool Anexa(const TCor& cor) = 0;
    virtual void Flush() = 0;
};

// ------------------------------------------------------------------------------------------------

class IArquivoSaida : public IDispositivoSaida
{
public:
    IArquivoSaida() = default;
    virtual ~IArquivoSaida() = default;

    virtual bool Aberto() const = 0;
    virtual const char* Caminho() const = 0;
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

class TWin32Viewport : public IDispositivoSaida
{
public:
    TWin32Viewport() = delete;
    TWin32Viewport(
        HDC hdc,
        std::int32_t x,
        std::int32_t y,
        std::int32_t w,
        std::int32_t h
    ) :
        hdc(hdc), xBeg(x), yBeg(y), w(w), h(h)
    {
    }

    bool Anexa(const TCor& cor) override
    {
        SetPixel(hdc, xBeg + x, yBeg + y, RGB(cor.R(), cor.G(), cor.B()));

        x++;
        if (x == w)
        {
            x = 0;
            y++;
        }

        return true;
    }

    void Flush() override
    {
        // NoOp
    }

private:
    HDC hdc;
    std::int32_t xBeg = 0;
    std::int32_t yBeg = 0;
    std::int32_t w = 0;
    std::int32_t h = 0;

    std::int32_t x = 0;
    std::int32_t y = 0;
};

// ------------------------------------------------------------------------------------------------

namespace FuncoesGerais
{
    bool EstaEm(double val, double min, double max, bool inicioFechado = true, bool fimFechado = true)
    {
        bool maiorQueMin = inicioFechado ? val >= min : val > min;
        bool menorQueMax = fimFechado ? val <= max : val < max;

        return maiorQueMin && menorQueMax;
    }

    uint8_t Trunca(double n)
    {
        return static_cast<uint8_t>(floor(n));
    }

    TCor Vec2Cor(const TVetor3D& v)
    {
        return { Trunca(v.X()), Trunca(v.Y()), Trunca(v.Z()) };
    }
    TVetor3D Cor2Vec(const TCor& c)
    {
        return { c.R() / 255.0, c.G() / 255.0, c.B() / 255.0 };
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

using TImagem = TMatriz<TCor, uint16_t>;

// ------------------------------------------------------------------------------------------------

class TTextura
{
public:
    TTextura()
    {
        _pixels = new TImagem(0, 0);
    }

    TTextura(const TTextura& outra)
    {
        _k = outra._k;
        _pixels = new TImagem(*outra._pixels);
    }

    TTextura(TTextura&& outra)
    {
        _k = outra._k;
        _pixels = outra._pixels;
        outra._pixels = nullptr;
    }

    TTextura& operator=(const TTextura& outra)
    {
        if (&outra != this)
        {
            _k = outra._k;
            _pixels = new TImagem(*outra._pixels);
        }

        return *this;
    }

    TTextura& operator=(TTextura& outra)
    {
        if (&outra != this)
        {
            _k = outra._k;
            _pixels = outra._pixels;
            outra._pixels = nullptr;
        }

        return *this;
    }

    virtual ~TTextura()
    {
        delete _pixels;
    }

    void Carrega(const std::string& caminho)
    {
        bmp::Bitmap img;
        img.load(caminho);

        const uint16_t w = static_cast<uint16_t>(img.width());
        const uint16_t h = static_cast<uint16_t>(img.height());

        delete _pixels;
        _pixels = new TImagem(h, w);

        for (std::size_t x = 0; x < w; x++)
        {
            for (std::size_t y = 0; y < h; y++)
            {
                const bmp::Pixel& pixel = img.get(x, y);
                (*_pixels)[y + 1][x + 1] = { pixel.r, pixel.g, pixel.b };
            }
        }
    }

    TCor Pixel(uint16_t x, uint16_t y) const
    {
        return (*_pixels)[y + 1][x + 1];
    }

    uint16_t Largura() const
    {
        return _pixels->NumeroColunas();
    }
    uint16_t Altura() const
    {
        return _pixels->NumeroLinhas();
    }

    double K() const
    {
        return _k;
    }
    void K(double k)
    {
        _k = k;
    }

private:
    double _k = 1.0;
    TImagem* _pixels = nullptr;
};

// ------------------------------------------------------------------------------------------------

class TMaterial
{
public:
    TMaterial() = default;

    TMaterial(const TMaterial& outro)
    {
        CopiaDadosPrimitivos(outro);

        if (outro._tex != nullptr)
        {
            _tex = new TTextura(*outro._tex);
        }
    }

    TMaterial(TMaterial&& outro)
    {
        CopiaDadosPrimitivos(outro);

        _tex = outro._tex;
        outro._tex = nullptr;
    }

    TMaterial& operator=(const TMaterial& outro)
    {
        if (&outro != this)
        {
            CopiaDadosPrimitivos(outro);

            if (outro._tex != nullptr)
            {
                _tex = new TTextura(*outro._tex);
            }
        }

        return *this;
    }

    TMaterial& operator=(TMaterial&& outro)
    {
        if (&outro != this)
        {
            CopiaDadosPrimitivos(outro);

            _tex = outro._tex;
            outro._tex = nullptr;
        }

        return *this;
    }

    virtual ~TMaterial()
    {
        delete _tex;
    }

    double M() const
    {
        return _m;
    }
    void M(double m)
    {
        _m = m;
    }

    double KdR(double u, double v) const
    {
        return _tex != nullptr ? KdTex(u, v).X() : _kdR;
    }
    double KdR() const
    {
        return _kdR;
    }
    void KdR(double r)
    {
        _kdR = r;
    }
    double KdG(double u, double v) const
    {
        return _tex != nullptr ? KdTex(u, v).Y() : _kdG;
    }
    double KdG() const
    {
        return _kdG;
    }
    void KdG(double g)
    {
        _kdG = g;
    }
    double KdB(double u, double v) const
    {
        return _tex != nullptr ? KdTex(u, v).Z() : _kdB;
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

    double KaR() const
    {
        return _kaR;
    }
    void KaR(double r)
    {
        _kaR = r;
    }
    double KaG() const
    {
        return _kaG;
    }
    void KaG(double g)
    {
        _kaG = g;
    }
    double KaB() const
    {
        return _kaB;
    }
    void KaB(double b)
    {
        _kaB = b;
    }

    const TTextura* Textura() const
    {
        return _tex;
    }
    TTextura* Textura()
    {
        return _tex;
    }
    void CarregaTextura(const std::string& caminho)
    {
        _tex = new TTextura;
        _tex->Carrega(caminho);
    }

private:
    void CopiaDadosPrimitivos(const TMaterial& outro)
    {
        _m = outro._m;

        _kdR = outro._kdR;
        _kdG = outro._kdG;
        _kdB = outro._kdB;

        _keR = outro._keR;
        _keG = outro._keG;
        _keB = outro._keB;

        _kaR = outro._kaR;
        _kaG = outro._kaG;
        _kaB = outro._kaB;
    }

    TVetor3D KdTex(double u, double v) const
    {
        if (_tex != nullptr)
        {
            auto x = static_cast<uint16_t>(u * (_tex->Largura() - 1));
            auto y = static_cast<uint16_t>((1.0 - v) * (_tex->Altura() - 1));

            return FuncoesGerais::Cor2Vec(_tex->Pixel(x, y));
        }

        return { 0.0, 0.0, 0.0 };
    }

    double _m = 0.0;

    double _kdR = 0.0;
    double _kdG = 0.0;
    double _kdB = 0.0;

    double _keR = 0.0;
    double _keG = 0.0;
    double _keB = 0.0;

    double _kaR = 0.0;
    double _kaG = 0.0;
    double _kaB = 0.0;

    TTextura* _tex = nullptr;
};

// ------------------------------------------------------------------------------------------------

class IFonteLuminosa
{
public:
    IFonteLuminosa() = default;
    virtual ~IFonteLuminosa() = default;

    virtual IFonteLuminosa* Copia() const = 0;
};

// ------------------------------------------------------------------------------------------------

class TFontePontual : public IFonteLuminosa
{
public:
    TFontePontual() = default;
    TFontePontual(const TPonto3D& posicao, const TVetor3D& intensidade)
        : _p(posicao), _i(intensidade) {}

    IFonteLuminosa* Copia() const override
    {
        return new TFontePontual(*this);
    }

    const TPonto3D& Posicao() const
    {
        return _p;
    }
    const TVetor3D& Intensidade() const
    {
        return _i;
    }

private:
    TPonto3D _p;
    TVetor3D _i;
};

// ------------------------------------------------------------------------------------------------

class IEntidade3D
{
public:
    IEntidade3D() = default;
    virtual ~IEntidade3D() = default;

    virtual IEntidade3D* Copia() const = 0;
    virtual std::string Rotulo() const = 0;
    virtual TMaterial Material(const TRaio3D& raio) const = 0;
    virtual std::pair<double, double> CoordenadasUV(const TPonto3D& p) const = 0;

    virtual TVetor3D Normal(const TPonto3D& p, const TRaio3D& raio) const = 0;
    virtual std::vector<double> Intersecoes(const TRaio3D& raio) const = 0;
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

    TVetor3D Normal(const TPonto3D& p1, const TPonto3D& p2, const TPonto3D& p3)
    {
        const TVetor3D r1 = p2 - p1;
        const TVetor3D r2 = p3 - p1;

        return r1.Vetorial(r2);
    }

    double ProdutoMisto(const TVetor3D& v1, const TVetor3D& v2, const TVetor3D& v3)
    {
        return v1.Dot(v2.Vetorial(v3));
    }

    std::vector<double> IntersecoesValidas(
        const IEntidade3D& entidade,
        const TRaio3D& raio
    )
    {
        std::vector<double> intersecoesValidas;

        auto TemIntersecao = [&intersecoesValidas](double intersecao)
        {
            for (double intersecaoValida : intersecoesValidas)
            {
                if (intersecao == intersecaoValida)
                {
                    return true;
                }
            }

            return false;
        };

        std::vector<double> intersecoes = entidade.Intersecoes(raio);
        for (double intersecao : intersecoes)
        {
            if (intersecao >= 0.0 && !TemIntersecao(intersecao))
            {
                intersecoesValidas.push_back(intersecao);
            }
        }

        return intersecoesValidas;
    }
}

// ------------------------------------------------------------------------------------------------

class TEntidadeComposta : public IEntidade3D
{
public:
    TEntidadeComposta() = default;

    TEntidadeComposta(const TEntidadeComposta& outra)
    {
        for (const std::unique_ptr<IEntidade3D>& entidade : outra._entidades)
        {
            _entidades.push_back(std::unique_ptr<IEntidade3D>(entidade->Copia()));
        }
    }

    IEntidade3D* Copia() const override
    {
        return new TEntidadeComposta(*this);
    }

    std::string Rotulo() const override
    {
        return _rotulo;
    }
    void Rotulo(const std::string& rotulo)
    {
        _rotulo = rotulo;
    }

    TMaterial Material(const TRaio3D& raio) const override
    {
        TMaterial material;

        const IEntidade3D* entidade = EntidadeInterceptada(raio);
        if (entidade != nullptr)
        {
            material = entidade->Material(raio);
        }

        return material;
    }

    std::pair<double, double> CoordenadasUV(const TPonto3D& p) const override
    {
        return std::pair<double, double>(0.0, 0.0);
    }

    TVetor3D Normal(const TPonto3D& p, const TRaio3D& raio) const override
    {
        TVetor3D normal;

        const IEntidade3D* entidade = EntidadeInterceptada(raio);
        if (entidade != nullptr)
        {
            normal = entidade->Normal(p, raio);
        }

        return normal;
    }
    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        const IEntidade3D* entidade = EntidadeInterceptada(raio);
        if (entidade != nullptr)
        {
            intersecoes = entidade->Intersecoes(raio);
        }

        return intersecoes;
    }

    void Insere(const IEntidade3D& entidade)
    {
        _entidades.push_back(std::unique_ptr<IEntidade3D>(entidade.Copia()));
    }

protected:
    IEntidade3D* EntidadeInterceptada(const TRaio3D& raio) const
    {
        double intersecaoMaisProximaObservador = std::numeric_limits<double>::max();
        IEntidade3D* entidadeMaisProximaObservador = nullptr;

        for (const std::unique_ptr<IEntidade3D>& entidade : _entidades)
        {
            const std::vector<double> intersecoesValidas = FuncoesGeometricas::IntersecoesValidas(*entidade, raio);
            if (!intersecoesValidas.empty())
            {
                const double intersecaoEntidadeMaisProximaObservador = intersecoesValidas[0];
                if (intersecaoEntidadeMaisProximaObservador < intersecaoMaisProximaObservador)
                {
                    intersecaoMaisProximaObservador = intersecaoEntidadeMaisProximaObservador;
                    entidadeMaisProximaObservador = entidade.get();
                }
            }
        }

        return entidadeMaisProximaObservador;
    }

    std::string _rotulo;

    std::vector<std::unique_ptr<IEntidade3D>> _entidades;
};

// ------------------------------------------------------------------------------------------------

class TPlano : public IEntidade3D
{
public:
    TPlano() = delete;
    TPlano(const TPonto3D& posicao, const TVetor3D& direcao) : _p(posicao), _n(direcao) {}

    IEntidade3D* Copia() const override
    {
        return new TPlano(*this);
    }

    std::string Rotulo() const override
    {
        return _rotulo;
    }
    void Rotulo(const std::string& rotulo)
    {
        _rotulo = rotulo;
    }

    TMaterial Material(const TRaio3D&) const override
    {
        return _material;
    }
    void Material(const TMaterial& material)
    {
        _material = material;
    }

    std::pair<double, double> CoordenadasUV(const TPonto3D& p) const override
    {
        double u = 0.0;
        double v = 0.0;

        if (const TTextura* tex = _material.Textura())
        {
            u = p.X() / tex->K();
            u = u - floor(u);

            v = p.Z() / tex->K();
            v = v - floor(v);
        }

        return std::pair<double, double>(u, v);
    }

    const TPonto3D& PontoReferencia() const
    {
        return _p;
    }
    const TVetor3D& Normal() const
    {
        return _n;
    }

    TVetor3D Normal(const TPonto3D&, const TRaio3D&) const override
    {
        return Normal();
    }

    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        const double dn = raio.Direcao().Dot(_n);
        if (dn != 0.0)
        {
            intersecoes.push_back(-(_n.Dot(raio.Origem() - _p) / dn));
        }

        return intersecoes;
    }

protected:
    std::string _rotulo;
    TMaterial _material;

    TPonto3D _p;
    TVetor3D _n;
};

// ------------------------------------------------------------------------------------------------

class TSuperficieCircular : public TPlano
{
public:
    TSuperficieCircular() = delete;
    TSuperficieCircular(const TPonto3D& posicao, const TVetor3D& direcao, double raio)
        : TPlano { posicao, direcao }, _r(raio) {}

    IEntidade3D* Copia() const override
    {
        return new TSuperficieCircular(*this);
    }

    const TPonto3D& PontoReferencia() const
    {
        return _p;
    }

    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        std::vector<double> intersecoesPlano = TPlano::Intersecoes(raio);
        for (double intersecaoPlano : intersecoesPlano)
        {
            if (TVetor3D { raio.Ponto(intersecaoPlano) - _p }.Norma() <= _r)
            {
                intersecoes.push_back(intersecaoPlano);
            }
        }

        return intersecoes;
    }

private:
    double _r;
};

// ------------------------------------------------------------------------------------------------

class TSuperficieTriangular : public TPlano
{
public:
    TSuperficieTriangular() = delete;
    TSuperficieTriangular(const TPonto3D& p1, const TPonto3D& p2, const TPonto3D& p3)
        : _p1(p1),
          _p2(p2),
          _p3(p3),
          TPlano { p1, FuncoesGeometricas::Normal(p1, p2, p3).Normalizado() }
    {
        _N = FuncoesGeometricas::Normal(p1, p2, p3);
    }

    IEntidade3D* Copia() const override
    {
        return new TSuperficieTriangular(*this);
    }

    const TPonto3D& PontoReferencia() const
    {
        return _p;
    }

    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        std::vector<double> intersecoesPlano = TPlano::Intersecoes(raio);
        for (double intersecaoPlano : intersecoesPlano)
        {
            const TPonto3D pI = raio.Ponto(intersecaoPlano);

            const TVetor3D s1 = _p1 - pI;
            const TVetor3D s2 = _p2 - pI;
            const TVetor3D s3 = _p3 - pI;

            const double normaN = _N.Norma();
            const double c1 = FuncoesGeometricas::ProdutoMisto(_n, s3, s1) / normaN;
            const double c2 = FuncoesGeometricas::ProdutoMisto(_n, s1, s2) / normaN;
            const double c3 = 1.0 - c1 - c2;
            const bool intersecaoDentro = c1 > 0.0 && c2 > 0.0 && c3 > 0.0;

            if (intersecaoDentro)
            {
                intersecoes.push_back(intersecaoPlano);
            }
        }

        return intersecoes;
    }

private:
    TPonto3D _p1;
    TPonto3D _p2;
    TPonto3D _p3;

    TVetor3D _N;
};

// ------------------------------------------------------------------------------------------------

class TSuperficieRetangular : public TEntidadeComposta
{
public:
    TSuperficieRetangular(
        const TPonto3D& p1,
        const TPonto3D& p2,
        const TPonto3D& p3,
        const TPonto3D& p4
    )
    {
        Insere(TSuperficieTriangular { p1, p2, p3 });
        Insere(TSuperficieTriangular { p1, p3, p4 });

        AtualizaReferencias();
    }

    TSuperficieRetangular(const TSuperficieRetangular& outra) : TEntidadeComposta(outra)
    {
        AtualizaReferencias();
    }

    IEntidade3D* Copia() const override
    {
        return new TSuperficieRetangular(*this);
    }

    void Material(const TMaterial& material)
    {
        _t1->Material(material);
        _t2->Material(material);
    }

private:
    void AtualizaReferencias()
    {
        _t1 = static_cast<TSuperficieTriangular*>(_entidades[0].get());
        _t2 = static_cast<TSuperficieTriangular*>(_entidades[1].get());
    }

    TSuperficieTriangular* _t1 = nullptr;
    TSuperficieTriangular* _t2 = nullptr;
};

// ------------------------------------------------------------------------------------------------

class TCubo : public TEntidadeComposta
{
public:
    TCubo(const TPonto3D& pCentro, double aresta) : _a(aresta), _sa(0.5 * aresta)
    {
        const TPonto3D pFaceInf1 { pCentro.X() - _sa, pCentro.Y(), pCentro.Z() + _sa };
        const TPonto3D pFaceInf2 { pCentro.X() + _sa, pCentro.Y(), pCentro.Z() + _sa };
        const TPonto3D pFaceInf3 { pCentro.X() + _sa, pCentro.Y(), pCentro.Z() - _sa };
        const TPonto3D pFaceInf4 { pCentro.X() - _sa, pCentro.Y(), pCentro.Z() - _sa };
        
        const TPonto3D pFaceSup1 { pFaceInf1.X(), pFaceInf1.Y() + _a, pFaceInf1.Z() };
        const TPonto3D pFaceSup2 { pFaceInf2.X(), pFaceInf2.Y() + _a, pFaceInf2.Z() };
        const TPonto3D pFaceSup3 { pFaceInf3.X(), pFaceInf3.Y() + _a, pFaceInf3.Z() };
        const TPonto3D pFaceSup4 { pFaceInf4.X(), pFaceInf4.Y() + _a, pFaceInf4.Z() };

        const TSuperficieRetangular faceInferior { pFaceInf1, pFaceInf2, pFaceInf3, pFaceInf4 };
        const TSuperficieRetangular faceSuperior { pFaceSup1, pFaceSup2, pFaceSup3, pFaceSup4 };
        const TSuperficieRetangular faceLateralEsquerda { pFaceInf1, pFaceSup1, pFaceSup4, pFaceInf4 };
        const TSuperficieRetangular faceLateralDireita { pFaceInf2, pFaceInf3, pFaceSup3, pFaceSup2 };
        const TSuperficieRetangular faceTraseira { pFaceInf3, pFaceInf4, pFaceSup4, pFaceSup3 };
        const TSuperficieRetangular faceFrontal { pFaceInf1, pFaceInf2, pFaceSup2, pFaceSup1 };

        Insere(faceInferior);
        Insere(faceSuperior);
        Insere(faceLateralEsquerda);
        Insere(faceLateralDireita);
        Insere(faceTraseira);
        Insere(faceFrontal);

        AtualizaReferencias();
    }

    TCubo(const TCubo& outro) : TEntidadeComposta(outro)
    {
        AtualizaReferencias();
    }

    IEntidade3D* Copia() const override
    {
        return new TCubo(*this);
    }

    void Material(const TMaterial& material)
    {
        _fInferior->Material(material);
        _fSuperior->Material(material);
        _fLateralEsq->Material(material);
        _fLateralDir->Material(material);
        _fTraseira->Material(material);
        _fFrontal->Material(material);
    }

private:
    void AtualizaReferencias()
    {
        _fInferior = static_cast<TSuperficieRetangular*>(_entidades[0].get());
        _fSuperior = static_cast<TSuperficieRetangular*>(_entidades[1].get());
        _fLateralEsq = static_cast<TSuperficieRetangular*>(_entidades[2].get());
        _fLateralDir = static_cast<TSuperficieRetangular*>(_entidades[3].get());
        _fTraseira = static_cast<TSuperficieRetangular*>(_entidades[4].get());
        _fFrontal = static_cast<TSuperficieRetangular*>(_entidades[5].get());
    }

    double _a = 0.0;
    double _sa = 0.0;

    TSuperficieRetangular* _fInferior = nullptr;
    TSuperficieRetangular* _fSuperior = nullptr;
    TSuperficieRetangular* _fLateralEsq = nullptr;
    TSuperficieRetangular* _fLateralDir = nullptr;
    TSuperficieRetangular* _fTraseira = nullptr;
    TSuperficieRetangular* _fFrontal = nullptr;
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

    std::string Rotulo() const override
    {
        return _rotulo;
    }
    void Rotulo(const std::string& rotulo)
    {
        _rotulo = rotulo;
    }

    TMaterial Material(const TRaio3D&) const override
    {
        return _material;
    }
    void Material(const TMaterial& material)
    {
        _material = material;
    }

    std::pair<double, double> CoordenadasUV(const TPonto3D& p) const override
    {
        return std::pair<double, double>(0.0, 0.0);
    }

    const TPonto3D& Centro() const { return _centro; }
    double Raio() const { return _raio; }

    TVetor3D Normal(const TPonto3D& p, const TRaio3D&) const override
    {
        return TVetor3D { p - Centro() } / Raio();
    }

    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        const TVetor3D w = raio.Origem() - Centro();
        const TVetor3D d = raio.Direcao();

        const double a = d.Dot(d);
        const double b = 2.0 * (w.Dot(d));
        const double c = w.Dot(w) - Raio() * Raio();

        const double delta = b * b - 4.0 * a * c;
        if (delta < 0.0)
        {
            return intersecoes;
        }
        intersecoes.push_back((-b - sqrt(delta)) / 2.0 * a);
        intersecoes.push_back((-b + sqrt(delta)) / 2.0 * a);

        return intersecoes;
    }

private:
    std::string _rotulo;
    TMaterial _material;

    TPonto3D _centro;
    double _raio = 0.0;
};

// ------------------------------------------------------------------------------------------------

class TSuperficieCilindrica : public IEntidade3D
{
public:
    TSuperficieCilindrica(const TPonto3D& pCentroBase, double raio, double altura, const TVetor3D& direcao)
        : _c(pCentroBase), _r(raio), _h(altura), _d(direcao.Normalizado()) {}

    IEntidade3D* Copia() const override
    {
        return new TSuperficieCilindrica(*this);
    }

    std::string Rotulo() const override
    {
        return _rotulo;
    }
    void Rotulo(const std::string& rotulo)
    {
        _rotulo = rotulo;
    }

    TMaterial Material(const TRaio3D&) const override
    {
        return _material;
    }
    void Material(const TMaterial& material)
    {
        _material = material;
    }

    std::pair<double, double> CoordenadasUV(const TPonto3D& p) const override
    {
        return std::pair<double, double>(0.0, 0.0);
    }

    // nos dois metodos abaixo, as notas de aula estavam bem diferentes
    // das fotos da lousa. inicialmente tentei implementar pelas fotos,
    // mas ocorreram diversos problemas, entao fui pelas notas, que
    // pareceu funcionar com alguns ajustes. depois avaliar o que errei
    // ao fazer a versao da lousa (p.s: aparentemente era bem mais
    // pesada, pelo uso ostensivo de matrizes)
    TVetor3D Normal(const TPonto3D& p, const TRaio3D&) const override
    {
        const TVetor3D s = p - _c;
        const TVetor3D s_ = _d * s.Dot(_d);

        return TVetor3D { s - s_ }.Normalizado();
    }
    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        const TVetor3D& dc = _d;
        const TVetor3D& dr = raio.Direcao();

        const TVetor3D aux = raio.Origem() - _c;
        const TVetor3D v = aux - dc * aux.Dot(dc);
        const TVetor3D w = dr - dc * dr.Dot(dc);

        const double a = w.Dot(w);
        const double b = 2.0 * v.Dot(w);
        const double c = v.Dot(v) - _r * _r;

        // talvez fosse bom avaliar "a" como fabs(a) < 1e-12
        const double delta = b * b - 4.0 * a * c;
        if (a == 0.0 || delta < 0.0)
        {
            return intersecoes;
        }

        // passar para o solver generico depois
        const double t1 = (-b - sqrt(delta)) / (2.0 * a);
        const double t2 = (-b + sqrt(delta)) / (2.0 * a);
        const TPonto3D p1 = raio.Ponto(t1);
        const TPonto3D p2 = raio.Ponto(t2);
        const TVetor3D s1 = p1 - _c;
        const TVetor3D s2 = p2 - _c;
        const double alturaT1 = s1.Dot(dc);
        const double alturaT2 = s2.Dot(dc);
        const bool t1Valido = FuncoesGerais::EstaEm(alturaT1, 0.0, _h);
        const bool t2Valido = FuncoesGerais::EstaEm(alturaT2, 0.0, _h);

        if (t1Valido)
        {
            intersecoes.push_back(t1);
        }
        if (t2Valido)
        {
            intersecoes.push_back(t2);
        }

        return intersecoes;
    }

    const TPonto3D& CentroBase() const
    {
        return _c;
    }
    double Raio() const
    {
        return _r;
    }
    double Altura() const
    {
        return _h;
    }
    const TVetor3D& Direcao() const
    {
        return _d;
    }

private:
    std::string _rotulo;
    TMaterial _material;

    TPonto3D _c;
    double _r;
    double _h;
    TVetor3D _d;
};

// ------------------------------------------------------------------------------------------------

class TCilindro : public TEntidadeComposta
{
public:
    TCilindro(const TPonto3D& pCentroBase, double raio, double altura, const TVetor3D& direcao)
    {
        Insere(TSuperficieCilindrica { pCentroBase, raio, altura, direcao });
        Insere(TSuperficieCircular { pCentroBase, direcao, raio });
        Insere(TSuperficieCircular { pCentroBase + direcao * altura, direcao, raio });

        _lateral = static_cast<TSuperficieCilindrica*>(_entidades[0].get());
        _base = static_cast<TSuperficieCircular*>(_entidades[1].get());
        _topo = static_cast<TSuperficieCircular*>(_entidades[2].get());
    }

    void Material(const TMaterial& material)
    {
        _lateral->Material(material);
        _base->Material(material);
        _topo->Material(material);
    }

    const TPonto3D& CentroBase() const
    {
        return _lateral->CentroBase();
    }
    double Raio() const
    {
        return _lateral->Raio();
    }
    double Altura() const
    {
        return _lateral->Altura();
    }
    const TVetor3D& Direcao() const
    {
        return _lateral->Direcao();
    }
private:
    TSuperficieCilindrica* _lateral = nullptr;
    TSuperficieCircular* _base = nullptr;
    TSuperficieCircular* _topo = nullptr;
};

// ------------------------------------------------------------------------------------------------

class TSuperficieConica : public IEntidade3D
{
public:
    TSuperficieConica(const TPonto3D& pCentroBase, double raioBase, double altura, const TVetor3D& direcao)
        : _c(pCentroBase), _r(raioBase), _h(altura), _d(direcao.Normalizado())
    {
        _v = _c + _d * _h;
        _cosTetaPow2 = (_h * _h) / (_h * _h + _r * _r);
    }

    IEntidade3D* Copia() const override
    {
        return new TSuperficieConica(*this);
    }

    std::string Rotulo() const override
    {
        return _rotulo;
    }
    void Rotulo(const std::string& rotulo)
    {
        _rotulo = rotulo;
    }

    TMaterial Material(const TRaio3D&) const override
    {
        return _material;
    }
    void Material(const TMaterial& material)
    {
        _material = material;
    }

    std::pair<double, double> CoordenadasUV(const TPonto3D& p) const override
    {
        return std::pair<double, double>(0.0, 0.0);
    }

    // nos dois metodos abaixo, as notas de aula estavam bem diferentes
    // das fotos da lousa. inicialmente tentei implementar pelas fotos,
    // mas ocorreram diversos problemas, entao fui pelas notas, que
    // pareceu funcionar com alguns ajustes. depois avaliar o que errei
    // ao fazer a versao da lousa (p.s: aparentemente era bem mais
    // pesada, pelo uso ostensivo de matrizes)
    TVetor3D Normal(const TPonto3D& p, const TRaio3D& r) const override
    {
        const TVetor3D w = p - _v;

        return TVetor3D { _d * _d.Dot(w) - w * _cosTetaPow2 }.Normalizado();
    }
    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        const TVetor3D& dr = raio.Direcao();
        const TVetor3D& dc = _d;

        const TVetor3D v = _v - raio.Origem();
        const double dn = dr.Dot(dc);
        const double dd = dr.Dot(dr);
        const double vn = v.Dot(dc);
        const double vd = v.Dot(dr);
        const double vv = v.Dot(v);

        const double a = dn * dn - dd * _cosTetaPow2;
        const double b = vd * _cosTetaPow2 - vn * dn;
        const double c = vn * vn - vv * _cosTetaPow2;

        // 1) ver aqui como escrever da forma tradicional, ex: 4ac,
        //    para poder passar para o solver generico depois
        // 2) talvez fosse bom avaliar "a" como fabs(a) < 1e-12
        const double delta = b * b - a * c;
        if (a == 0.0 || delta < 0.0)
        {
            return intersecoes;
        }

        // ver aqui como escrever da forma tradicional, ex: / 2.0a
        // para poder passar para o solver generico depois
        const double t1 = (-b - sqrt(delta)) / (a);
        const double t2 = (-b + sqrt(delta)) / (a);
        const TPonto3D p1 = raio.Ponto(t1);
        const TPonto3D p2 = raio.Ponto(t2);
        const double proj1 = dc.Dot(_v - p1);
        const double proj2 = dc.Dot(_v - p2);
        const bool t1Valido = FuncoesGerais::EstaEm(proj1, 0.0, _h);
        const bool t2Valido = FuncoesGerais::EstaEm(proj2, 0.0, _h);

        if (t1Valido)
        {
            intersecoes.push_back(t1);
        }
        if (t2Valido)
        {
            intersecoes.push_back(t2);
        }

        // nao entendi pq foi necessario ordenar, acho que ja devia vir
        // na ordem certa, mas pelo visto resolveu, entao devia estar
        // vindo um t maior no comeco, avaliar pq
        std::sort(intersecoes.begin(), intersecoes.end());

        return intersecoes;
    }

    const TPonto3D& CentroBase() const
    {
        return _c;
    }
    double RaioBase() const
    {
        return _r;
    }
    double Altura() const
    {
        return _h;
    }
    const TVetor3D& Direcao() const
    {
        return _d;
    }

private:
    std::string _rotulo;
    TMaterial _material;

    TPonto3D _c;
    double _r;
    double _h;
    TVetor3D _d;

    TPonto3D _v;
    double _cosTetaPow2 = 0.0;
};

// ------------------------------------------------------------------------------------------------

class TCone : public TEntidadeComposta
{
public:
    TCone(const TPonto3D& pCentroBase, double raioBase, double altura, const TVetor3D& direcao)
    {
        Insere(TSuperficieConica { pCentroBase, raioBase, altura, direcao });
        Insere(TSuperficieCircular { pCentroBase, direcao, raioBase });

        _lateral = static_cast<TSuperficieConica*>(_entidades[0].get());
        _base = static_cast<TSuperficieCircular*>(_entidades[1].get());
    }

    void Material(const TMaterial& material)
    {
        _lateral->Material(material);
        _base->Material(material);
    }

    const TPonto3D& CentroBase() const
    {
        return _lateral->CentroBase();
    }
    double RaioBase() const
    {
        return _lateral->RaioBase();
    }
    double Altura() const
    {
        return _lateral->Altura();
    }
    const TVetor3D& Direcao() const
    {
        return _lateral->Direcao();
    }
private:
    TSuperficieConica* _lateral = nullptr;
    TSuperficieCircular* _base = nullptr;
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

    double IambR() const
    {
        return _iAmbR;
    }
    void IambR(double iAmbR)
    {
        _iAmbR = iAmbR;
    }
    double IambG() const
    {
        return _iAmbG;
    }
    void IambG(double iAmbG)
    {
        _iAmbG = iAmbG;
    }
    double IambB() const
    {
        return _iAmbB;
    }
    void IambB(double iAmbB)
    {
        _iAmbB = iAmbB;
    }

    void Insere(const IEntidade3D& entidade)
    {
        _entidades.push_back(std::unique_ptr<IEntidade3D>(entidade.Copia()));
    }
    void Insere(const IFonteLuminosa& fonte)
    {
        _fontes.push_back(std::unique_ptr<IFonteLuminosa>(fonte.Copia()));
    }

    void Renderizar(IDispositivoSaida& arq)
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
        const TVetor3D d = FuncoesGeometricas::Versor(_p0, p);
        const TRaio3D raio { _p0, d };

        double intersecaoMaisProximaObservador = std::numeric_limits<double>::max();
        IEntidade3D* entidadeMaisProximaObservador = nullptr;

        for (const std::unique_ptr<IEntidade3D>& entidade : _entidades)
        {
            const std::vector<double> intersecoesValidas = FuncoesGeometricas::IntersecoesValidas(*entidade, raio);
            if (!intersecoesValidas.empty())
            {
                const double intersecaoEntidadeMaisProximaObservador = intersecoesValidas[0];
                if (intersecaoEntidadeMaisProximaObservador < intersecaoMaisProximaObservador)
                {
                    intersecaoMaisProximaObservador = intersecaoEntidadeMaisProximaObservador;
                    entidadeMaisProximaObservador = entidade.get();
                }
            }
        }

        TCor pixel = _bgColor;

        if (entidadeMaisProximaObservador != nullptr)
        {
            pixel = Cor(*entidadeMaisProximaObservador, raio, intersecaoMaisProximaObservador);
        }

        return pixel;
    }

    TCor Cor(const IEntidade3D& entidade, const TRaio3D& raio, double ti) const
    {   
        const TPonto3D pi = raio.Ponto(ti);
        const std::pair<double, double> uv = entidade.CoordenadasUV(pi);

        const TMaterial& material = entidade.Material(raio);
        const double kdR = material.KdR(uv.first, uv.second);
        const double kdG = material.KdG(uv.first, uv.second);
        const double kdB = material.KdB(uv.first, uv.second);
        const TPonto3D kd { kdR, kdG, kdB };
        const TPonto3D ke { material.KeR(), material.KeG(), material.KeB() };
        const TPonto3D ka { material.KaR(), material.KaG(), material.KaB() };

        const TVetor3D iAmb { IambR(), IambG(), IambB() };
        const auto ia = iAmb.Arroba(ka);

        const TVetor3D n = entidade.Normal(pi, raio);
        const TVetor3D v = raio.Direcao() * -1.0;

        TVetor3D i = ia;
        for (const std::unique_ptr<IFonteLuminosa>& fonte : _fontes)
        {
            if (auto fontePontual = dynamic_cast<const TFontePontual*>(fonte.get()))
            {
                const TPonto3D& pFonte = fontePontual->Posicao();
                const TVetor3D& iFonte = fontePontual->Intensidade();
                
                const TVetor3D L = pFonte - pi;
                const TVetor3D l = L.Normalizado();
                const TRaio3D shadowRay { pi, l };

                bool temEntidadeBloqueandoLuz = false;
                for (const std::unique_ptr<IEntidade3D>& outraEntidade : _entidades)
                {
                    if (&entidade != outraEntidade.get())
                    {
                        const std::vector<double> intersecoes = FuncoesGeometricas::IntersecoesValidas(
                            *outraEntidade, shadowRay
                        );
                        
                        if (!intersecoes.empty())
                        {
                            const double intersecaoMaisProxima = intersecoes[0];

                            if (intersecaoMaisProxima < L.Norma())
                            {
                                temEntidadeBloqueandoLuz = true;
                            }
                        }
                    }
                }

                if (!temEntidadeBloqueandoLuz)
                {
                    const TVetor3D r = (n * 2.0 * n.Dot(l)) - l;
                    const double fd = std::max(0.0, n.Dot(l));
                    const double fe = pow(std::max(v.Dot(r), 0.0), material.M());

                    const auto id = iFonte.Arroba(kd) * fd;
                    const auto ie = iFonte.Arroba(ke) * fe;

                    const TVetor3D iCorrente = id + ie;
                    i += iCorrente;
                }
            }
        }

        TVetor3D c = i.Arroba(TVetor3D(255.0, 255.0, 255.0));
        c.Clamp(0.0, 255.0);

        return FuncoesGerais::Vec2Cor(c);
    }

    TArquivoLOG* _arqLog = nullptr;

    TJanela _janela;
    TPonto3D _p0; // olho do pintor (origem)

    TCor _bgColor;
    double _iAmbR = 0.0;
    double _iAmbG = 0.0;
    double _iAmbB = 0.0;

    std::vector<std::unique_ptr<IEntidade3D>> _entidades;
    std::vector<std::unique_ptr<IFonteLuminosa>> _fontes;
};

// ------------------------------------------------------------------------------------------------

TMaterial FabricaMaterialHomogeneo(const TVetor3D& k, double m)
{
    TMaterial material;

    material.KdR(k.X());
    material.KdG(k.Y());
    material.KdB(k.Z());
    material.KeR(k.X());
    material.KeG(k.Y());
    material.KeB(k.Z());
    material.KaR(k.X());
    material.KaG(k.Y());
    material.KaB(k.Z());
    material.M(m);

    return material;
}

// ------------------------------------------------------------------------------------------------

TPlano FabricaChao()
{
    TMaterial material;
    material.KdR(0.0);
    material.KdG(0.0);
    material.KdB(0.0);
    material.KeR(0.0);
    material.KeG(0.0);
    material.KeB(0.0);
    material.KaR(0.0);
    material.KaG(0.0);
    material.KaB(0.0);
    material.M(1.0);
    material.CarregaTextura("tex.bmp");
    material.Textura()->K(30.0);

    TPlano planoChao { { 0.0, -150.0, 0.0 }, { 0.0, 1.0, 0.0 } };
    planoChao.Rotulo("PLANO_CHAO");
    planoChao.Material(material);

    return planoChao;
}

// ------------------------------------------------------------------------------------------------

TPlano FabricaParedeLateralDireita()
{
    TPlano planoParedeDireita { { 200.0, -150.0, 0.0 }, { -1.0, 0.0, 0.0 } };
    planoParedeDireita.Rotulo("PLANO_PAREDE_LATERAL_DIR");
    planoParedeDireita.Material(FabricaMaterialHomogeneo({ 0.686, 0.933, 0.933 }, 1.0));

    return planoParedeDireita;
}

// ------------------------------------------------------------------------------------------------

TPlano FabricaParedeFrontal()
{
    TPlano planoParedeFrontal { { 200.0, -150.0, -400.0 }, { 0.0, 0.0, 1.0 } };
    planoParedeFrontal.Rotulo("PLANO_PAREDE_FRONTAL");
    planoParedeFrontal.Material(FabricaMaterialHomogeneo({ 0.686, 0.933, 0.933 }, 1.0));

    return planoParedeFrontal;
}

// ------------------------------------------------------------------------------------------------

TPlano FabricaParedeLateralEsquerda()
{
    TPlano planoParedeEsquerda { { -200.0, -150.0, 0.0 }, { 1.0, 0.0, 0.0 } };
    planoParedeEsquerda.Rotulo("PLANO_PAREDE_LATERAL_ESQ");
    planoParedeEsquerda.Material(FabricaMaterialHomogeneo({ 0.686, 0.933, 0.933 }, 1.0));

    return planoParedeEsquerda;
}

// ------------------------------------------------------------------------------------------------

TPlano FabricaTeto()
{
    TPlano planoTeto { { 0.0, 150.0, 0.0 }, { 0.0, -1.0, 0.0 } };
    planoTeto.Rotulo("PLANO_TETO");
    planoTeto.Material(FabricaMaterialHomogeneo({ 0.933, 0.933, 0.933 }, 1.0));

    return planoTeto;
}

// ------------------------------------------------------------------------------------------------

TCilindro FabricaCilindro()
{
    const TPonto3D cBaseCilindro { 0.0, -150.0, -200.0 };
    const double rBaseCilindro = 5.0;
    const double hCilindro = 90.0;
    const TVetor3D dCilindro { 0.0, 1.0, 0.0 };

    TCilindro cilindro { cBaseCilindro, rBaseCilindro, hCilindro, dCilindro };
    cilindro.Rotulo("CILINDRO_1");
    cilindro.Material(FabricaMaterialHomogeneo({ 0.824, 0.706, 0.549 }, 10.0));

    return cilindro;
}

// ------------------------------------------------------------------------------------------------

TCone FabricaCone()
{
    const TPonto3D cBaseCone { 0.0, -60.0, -200.0 };
    const double rBaseCone = 90.0;
    const double hCone = 150.0;
    const TVetor3D& dCone = { 0.0, 1.0, 0.0 };

    TCone cone { cBaseCone, rBaseCone, hCone, dCone };
    cone.Rotulo("CONE_1");
    cone.Material(FabricaMaterialHomogeneo({ 0.0, 1.0, 0.498 }, 10.0));

    return cone;
}

// ------------------------------------------------------------------------------------------------

TCubo FabricaCubo()
{
    TCubo cubo { { 0.0, -150.0, -165.0 }, 40.0 };
    cubo.Rotulo("CUBO_1");
    cubo.Material(FabricaMaterialHomogeneo({ 1.0, 0.078, 0.576 }, 10.0));

    return cubo;
}

// ------------------------------------------------------------------------------------------------

TEsfera FabricaEsfera()
{
    TEsfera esfera { { 0.0, 95.0, -200.0 }, 5.0 };
    esfera.Rotulo("ESFERA_1");
    esfera.Material(FabricaMaterialHomogeneo({ 0.854, 0.647, 0.125 }, 10.0));

    return esfera;
}

// ------------------------------------------------------------------------------------------------

TFontePontual FabricaFontePontual()
{
    return { { -100.0, 140.0, -20.0 }, { 0.7, 0.7, 0.7 } };
}

// ------------------------------------------------------------------------------------------------

TCena3D FabricaCena()
{
    const TPonto3D p0 { 0.0, 0.0, 0.0 };

    const double wJanela = 60.0;
    const double hJanela = 60.0;
    const double dJanela = 30.0;
    const uint16_t wCanvas = 500u;
    const uint16_t hCanvas = 500u;
    const TJanela janela { { 0.0, 0.0, -dJanela }, wJanela, hJanela, wCanvas, hCanvas };

    TCena3D cena { p0, janela };
    cena.BgColor({ 100u, 100u, 100u });
    cena.IambR(0.3);
    cena.IambG(0.3);
    cena.IambB(0.3);

    cena.Insere(FabricaChao());
    cena.Insere(FabricaParedeLateralDireita());
    cena.Insere(FabricaParedeFrontal());
    cena.Insere(FabricaParedeLateralEsquerda());
    cena.Insere(FabricaTeto());
    cena.Insere(FabricaCilindro());
    cena.Insere(FabricaCone());
    cena.Insere(FabricaCubo());
    cena.Insere(FabricaEsfera());
    cena.Insere(FabricaFontePontual());

    return cena;
}

// ------------------------------------------------------------------------------------------------

std::unique_ptr<IArquivoSaida> FabricaArquivo(
    const std::string& nome,
    EFormatoImagem formato,
    const TCena3D& cena
)
{
    const TJanela janela = cena.Janela();

    return FuncoesGerais::FabricaArquivo(
        formato, nome, janela.LarguraCanvas(), janela.AlturaCanvas()
    );
}

// ------------------------------------------------------------------------------------------------

bool RenderizaImagem()
{
    TCena3D cena = FabricaCena();
    const std::unique_ptr<IArquivoSaida> arq = FabricaArquivo("a", EFormatoImagem::BMP, cena);

    const bool erro = !arq->Aberto();
    if (!erro)
    {
        cena.Renderizar(*arq);
    }

    return erro;
}

// ------------------------------------------------------------------------------------------------

class TMainWindow
{
public:
    TMainWindow(
        HINSTANCE hInstance,
        PWSTR pCmdLine,
        int nCmdShow
    ) :
        _hInstance(hInstance),
        _pCmdLine(pCmdLine),
        _nCmdShow(nCmdShow)
    {
    };

    TMainWindow& NomeClasse(const std::wstring& nomeClasse)
    {
        _clsName = nomeClasse;
        return *this;
    };

    TMainWindow& TituloJanela(const std::wstring& tituloJanela)
    {
        _wndTitle = tituloJanela;
        return *this;
    };

    void Executa()
    {
        if (MontaJanelaPrincipal())
        {
            ExibeJanelaPrincipal();
            ProcessaMensagens();
        }
    }

private:
    bool MontaJanelaPrincipal()
    {
        _wndCls = CriaClasse();
        RegisterClass(&_wndCls);

        _hWnd = CriaJanelaPrincipal();
        _hWndBtnRenderScene = CriaBotao(_hWnd, _IDC_BTN_RENDER_SCENE, L"Renderiza Cena", 16, 16, 128, 32);

        return _hWnd != nullptr && _hWndBtnRenderScene != nullptr;
    }

    WNDCLASS CriaClasse() const
    {
        WNDCLASS wc = {};
        wc.lpfnWndProc   = TMainWindow::WindowProc;
        wc.hInstance     = _hInstance;
        wc.lpszClassName = _clsName.c_str();

        return wc;
    }

    HWND CriaJanelaPrincipal() const
    {
        HWND hWnd = CreateWindowEx(
            0,                              // Optional window styles
            _clsName.c_str(),               // Window class
            _wndTitle.c_str(),              // Window text
            WS_OVERLAPPEDWINDOW,            // Window style

            // Size and position
            CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT,

            NULL,       // Parent window    
            NULL,       // Menu
            _hInstance, // Instance handle
            NULL        // Additional application data
        );

        return hWnd;
    }

    HWND CriaBotao(HWND parent, int id, const std::wstring& txt, int x, int y, int w, int h) const
    {
        HWND hWnd = CreateWindowEx( 
            0,                   // Optional window styles
            L"BUTTON",           // Predefined class; Unicode assumed
            txt.c_str(),         // Button text
            WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON, // Styles
            x,                   // x position
            y,                   // y position
            w,                   // Button width
            h,                   // Button height
            parent,              // Parent window
            (HMENU)(intptr_t)id, // Control ID for child windows
            (HINSTANCE)GetWindowLongPtr(parent, GWLP_HINSTANCE), 
            NULL);               // Pointer not needed.

        return hWnd;
    }

    void ExibeJanelaPrincipal()
    {
        ShowWindow(_hWnd, _nCmdShow);
    }

    static void EvBtnRenderScene(HWND hWndJanelaPrincipal)
    {
        rendering = true;

        InvalidateRect(hWndJanelaPrincipal, NULL, true);
        UpdateWindow(hWndJanelaPrincipal);
    }

    static void RenderScene(HDC hdc)
    {
        auto cena = FabricaCena();
        const std::uint16_t x = 160;
        const std::uint16_t y = 16;
        const std::uint16_t w = cena.Janela().LarguraCanvas();
        const std::uint16_t h = cena.Janela().AlturaCanvas();

        RECT r { x, y, x + w, y + h };
        FillRect(hdc, &r, (HBRUSH)GetStockObject(BLACK_BRUSH));

        if (rendering)
        {
            HCURSOR cursorAntigo = GetCursor();
            SetCursor(LoadCursor(NULL, IDC_WAIT));

            TWin32Viewport wndVp(hdc, x, y, w, h);
            cena.Renderizar(wndVp);

            SetCursor(cursorAntigo);
        }

        rendering = false;
    }

    void ProcessaMensagens()
    {
        MSG msg = {};
        while (GetMessage(&msg, NULL, 0, 0))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

    static LRESULT CALLBACK WindowProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
    {
        switch (uMsg)
        {
            case WM_DESTROY   : return EvQuit();
            case WM_PAINT     : return EvPaint(hWnd);
            case WM_SETCURSOR : if (EvSetCursor(lParam)) return true; break;
            case WM_MOUSEMOVE : EvMouseMove(hWnd, lParam); break;
            case WM_COMMAND   : return EvCommand(hWnd, wParam);
        }

        return DefWindowProc(hWnd, uMsg, wParam, lParam);
    }

    static bool EvQuit()
    {
        PostQuitMessage(0);

        return 0;
    }

    static bool EvPaint(HWND hWnd)
    {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hWnd, &ps);
        FillRect(hdc, &ps.rcPaint, (HBRUSH) (COLOR_WINDOW+1));
        RenderScene(hdc);
        EndPaint(hWnd, &ps);

        return 0;
    }

    static bool EvSetCursor(LPARAM lParam)
    {
        bool ok = false;
        if (LOWORD(lParam) == HTCLIENT)
        {
            SetCursor(LoadCursor(NULL, IDC_ARROW));
            ok = true;
        }

        return ok;
    }

    static bool EvMouseMove(HWND hWnd, LPARAM lParam)
    {
        const int xMousePos = GET_X_LPARAM(lParam);
        const int yMousePos = GET_Y_LPARAM(lParam);

        wchar_t buf[32];
        GetClassName(hWnd, buf, 32);

        std::wstringstream ss;
        ss << buf << " ( X: " << xMousePos << ", Y: " << yMousePos << " )";
        auto str = ss.str();
        
        SetWindowText(hWnd, str.c_str());

        return true;
    }

    static bool EvCommand(HWND hWnd, WPARAM wParam)
    {
        if (HIWORD(wParam) == BN_CLICKED)
        {
            if (LOWORD(wParam) == _IDC_BTN_RENDER_SCENE)
            {
                EvBtnRenderScene(hWnd);
            }
        }

        return 0;
    }

    std::wstring _clsName;
    std::wstring _wndTitle;
    HINSTANCE _hInstance;
    PWSTR _pCmdLine;
    int _nCmdShow;

    HWND _hWnd;
    WNDCLASS _wndCls;

    HWND _hWndBtnRenderScene; static constexpr int _IDC_BTN_RENDER_SCENE = 1;

    static bool rendering;
};

bool TMainWindow::rendering = false;

// ------------------------------------------------------------------------------------------------

int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE /*hPrevInstance*/, PWSTR pCmdLine, int nCmdShow)
{
    const std::wstring tituloJanela = L"CG1";
    TMainWindow(hInstance, pCmdLine, nCmdShow)
        .NomeClasse(tituloJanela)
        .TituloJanela(tituloJanela)
        .Executa();

    return 0;
}

// ------------------------------------------------------------------------------------------------
