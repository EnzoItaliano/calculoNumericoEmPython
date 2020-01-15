
# calculoNumericoEmPython

Este projeto está relacionado ao prof. Wellington José Corrêa que leciona cálculo numérico na UTFPR de Campo Mourão e consiste em traduzir os códigos que são atualmente utilizados no wxMaxima para python fazendo com que os alunos tenham contato com uma linguagem de programação que já viram antes na disciplina de programação.

## Instalação

Para a instalação do **python3** em seu computador basta seguir os passos do tutorial de acordo com seu sistema operacional.

-  [Linux](https://python.org.br/instalacao-linux/)

-  [Windows](https://python.org.br/instalacao-windows/)

-  [Mac](https://python.org.br/instalacao-mac/)

É preciso usar o módulo **pip** para instalar dependências de algumas funções. A instalação é da seguinte forma:

- Linux
Simplesmente digite o comando.
	```
	sudo apt install python3-pip
	```
	Para verificar se foi instalado corretamente pode ser utilizado os passos descritos no tópico do Windows logo abaixo.

- Windows
No windows o pip é instalado junto ao python3 no instalador.

- Mac
	Digite os seguintes comandos.
	```
	curl -O https://bootstrap.pypa.io/get-pip.py
	```
	```
	sudo python3 get-pip.py
	```
**Verificação**
Para verificar se foi instalado corretamente pode se utilizar dos seguintes passos.
1. Digite no terminal.
	```
	pip3 install --user pybin
	```
2. Se der tudo certo digite.
	```
	python3 -m pybin
	```
	Isto irá dizer se o caminho foi bem configurado. Se uma mensagem for mostrado parecida com...
	```
	The user bin directory `/home/cefn/.local/bin` is not in PATH
	```
3. Digite o seguinte.
	```
	python3 -m pybin put
	```	
	Isto vai fazer com que programas instalados com o pip fiquem disponíveis para serem usados.

**Dependências**
- Matplotlib
```
python3 -m pip install -U matplotlib
```
- PrettyTable
```
pip3 install PrettyTable
```
- [SymPy](https://docs.sympy.org/1.5.1/install.html)
- PrettyMatrix
```
pip3 install prettymatrix
```