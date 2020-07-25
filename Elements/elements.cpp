/*****************************************************************************
© 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and 
to permit others to do so.


This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
    
    1.  Redistributions of source code must retain the above copyright notice, this list of 
        conditions and the following disclaimer.
 
    2.  Redistributions in binary form must reproduce the above copyright notice, this list of 
        conditions and the following disclaimer in the documentation and/or other materials 
        provided with the distribution.
 
    3.  Neither the name of the copyright holder nor the names of its contributors may be used 
        to endorse or promote products derived from this software without specific prior 
        written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

**********************************************************************************************/



#include <iostream>  // std::cout etc.
#include <cmath>

#include "utilities.h"
#include "elements.h"



using namespace utils;


namespace elements{

//==============================================================================
//   Function Definitions
//==============================================================================

/*

==========================
Representative Local Cell 
==========================

              K
              ^         J
              |        /
              |       /
                     /
      6------------------7
     /|                 /|
    / |                / |
   /  |               /  |
  /   |              /   | 
 /    |             /    |
4------------------5     |
|     |            |     | ----> I
|     |            |     |  
|     |            |     |
|     |            |     |
|     2------------|-----3
|    /             |    /
|   /              |   /
|  /               |  /
| /                | /         
|/                 |/
0------------------1


face 0: [0,1,3,2]
face 1: [4,5,7,6]
face 2: [0,1,5,4]
face 3: [2,3,7,6]
face 4: [0,2,6,4]
face 6; [1,3,7,5]

*/

// creates nodal positions with Labatto spacing
void labatto_nodes_1D(
    c_array_t <real_t> &lab_nodes_1D,
    const int &num){
    if (num == 1){
        lab_nodes_1D(0) = 0.0;
    }
    else if (num == 2){
        lab_nodes_1D(0) = -1.0;
        lab_nodes_1D(1) =  1.0;
    }
    else if (num == 3){
        lab_nodes_1D(0) = -1.0;
        lab_nodes_1D(1) =  0.0;
        lab_nodes_1D(2) =  1.0;
    }
    else if (num == 4){
        lab_nodes_1D(0) = -1.0;
        lab_nodes_1D(1) = -1.0/5.0*sqrt(5.0);
        lab_nodes_1D(2) =  1.0/5.0*sqrt(5.0);
        lab_nodes_1D(3) =  1.0;
    }
    else if (num == 5){
        lab_nodes_1D(0) = -1.0;
        lab_nodes_1D(1) = -1.0/7.0*sqrt(21.0);
        lab_nodes_1D(2) =  0.0;
        lab_nodes_1D(3) =  1.0/7.0*sqrt(21.0);
        lab_nodes_1D(4) =  1.0;
    }
    else if (num == 6){
        lab_nodes_1D(0) = -1.0;
        lab_nodes_1D(1) = -sqrt(1.0/21.0*(7.0 + 2.0*sqrt(7.0)));
        lab_nodes_1D(2) = -sqrt(1.0/21.0*(7.0 - 2.0*sqrt(7.0)));
        lab_nodes_1D(3) =  sqrt(1.0/21.0*(7.0 - 2.0*sqrt(7.0)));
        lab_nodes_1D(4) =  sqrt(1.0/21.0*(7.0 +2.0*sqrt(7.0)));
        lab_nodes_1D(5) =  1.0;
    }
    else if (num == 7){
        lab_nodes_1D(0) =  - 1.0E+00;
        lab_nodes_1D(1) =  - 0.830223896278566929872032213967E+00;
        lab_nodes_1D(2) =  - 0.468848793470714213803771881909E+00;
        lab_nodes_1D(3) =    0.0E+00;
        lab_nodes_1D(4) =    0.468848793470714213803771881909E+00;
        lab_nodes_1D(5) =    0.830223896278566929872032213967E+00;
        lab_nodes_1D(6) =    1.0E+00;
    }
    else if (num == 8){
        lab_nodes_1D(0) =  - 1.0E+00;
        lab_nodes_1D(1) =  - 0.871740148509606615337445761221E+00;
        lab_nodes_1D(2) =  - 0.591700181433142302144510731398E+00;
        lab_nodes_1D(3) =  - 0.209299217902478868768657260345E+00;
        lab_nodes_1D(4) =    0.209299217902478868768657260345E+00;
        lab_nodes_1D(5) =    0.591700181433142302144510731398E+00;
        lab_nodes_1D(6) =    0.871740148509606615337445761221E+00;
        lab_nodes_1D(7) =    1.0E+00;
    }
    else if (num == 9){
        lab_nodes_1D(0) =  - 1.0E+00;
        lab_nodes_1D(1) =  - 0.899757995411460157312345244418E+00;
        lab_nodes_1D(2) =  - 0.677186279510737753445885427091E+00;
        lab_nodes_1D(3) =  - 0.363117463826178158710752068709E+00;
        lab_nodes_1D(4) =    0.0E+00;
        lab_nodes_1D(5) =    0.363117463826178158710752068709E+00;
        lab_nodes_1D(6) =    0.677186279510737753445885427091E+00;
        lab_nodes_1D(7) =    0.899757995411460157312345244418E+00;
        lab_nodes_1D(8) =    1.0E+00;
    }
    else if (num == 10){
        lab_nodes_1D(0) =  - 1.0E+00;
        lab_nodes_1D(1) =  - 0.919533908166458813828932660822E+00;
        lab_nodes_1D(2) =  - 0.738773865105505075003106174860E+00;
        lab_nodes_1D(3) =  - 0.477924949810444495661175092731E+00;
        lab_nodes_1D(4) =  - 0.165278957666387024626219765958E+00;
        lab_nodes_1D(5) =    0.165278957666387024626219765958E+00;
        lab_nodes_1D(6) =    0.477924949810444495661175092731E+00;
        lab_nodes_1D(7) =    0.738773865105505075003106174860E+00;
        lab_nodes_1D(8) =    0.919533908166458813828932660822E+00;
        lab_nodes_1D(9) =   1.0E+00;

    }
    else if (num == 11){
        lab_nodes_1D(0) =  - 1.0E+00;
        lab_nodes_1D(1) =  - 0.934001430408059134332274136099E+00;
        lab_nodes_1D(2) =  - 0.784483473663144418622417816108E+00;
        lab_nodes_1D(3) =  - 0.565235326996205006470963969478E+00;
        lab_nodes_1D(4) =  - 0.295758135586939391431911515559E+00;
        lab_nodes_1D(5) =    0.0E+00;
        lab_nodes_1D(6) =    0.295758135586939391431911515559E+00;
        lab_nodes_1D(7) =    0.565235326996205006470963969478E+00;
        lab_nodes_1D(8) =    0.784483473663144418622417816108E+00;
        lab_nodes_1D(9) =   0.934001430408059134332274136099E+00;
        lab_nodes_1D(10) =   1.0E+00;
    }

    else if (num == 12){
        lab_nodes_1D(0) =  - 1.0E+00;
        lab_nodes_1D(1) =  - 0.944899272222882223407580138303E+00;
        lab_nodes_1D(2) =  - 0.819279321644006678348641581717E+00;
        lab_nodes_1D(3) =  - 0.632876153031869677662404854444E+00;
        lab_nodes_1D(4) =  - 0.399530940965348932264349791567E+00;
        lab_nodes_1D(5) =  - 0.136552932854927554864061855740E+00;
        lab_nodes_1D(6) =    0.136552932854927554864061855740E+00;
        lab_nodes_1D(7) =    0.399530940965348932264349791567E+00;
        lab_nodes_1D(8) =    0.632876153031869677662404854444E+00;
        lab_nodes_1D(9) =   0.819279321644006678348641581717E+00;
        lab_nodes_1D(10) =   0.944899272222882223407580138303E+00;
        lab_nodes_1D(11) =   1.0E+00;
    }

    else if (num == 13){
        lab_nodes_1D(0) =  - 1.0E+00;
        lab_nodes_1D(1) =  - 0.953309846642163911896905464755E+00;
        lab_nodes_1D(2) =  - 0.846347564651872316865925607099E+00;
        lab_nodes_1D(3) =  - 0.686188469081757426072759039566E+00;
        lab_nodes_1D(4) =  - 0.482909821091336201746937233637E+00;
        lab_nodes_1D(5) =  - 0.249286930106239992568673700374E+00;
        lab_nodes_1D(6) =    0.0E+00;
        lab_nodes_1D(7) =    0.249286930106239992568673700374E+00;
        lab_nodes_1D(8) =    0.482909821091336201746937233637E+00;
        lab_nodes_1D(9) =   0.686188469081757426072759039566E+00;
        lab_nodes_1D(10) =   0.846347564651872316865925607099E+00;
        lab_nodes_1D(11) =   0.953309846642163911896905464755E+00;
        lab_nodes_1D(12) =   1.0E+00;
    }

    else if (num == 14){
        lab_nodes_1D(0) =  - 1.0E+00;
        lab_nodes_1D(1) =  - 0.959935045267260901355100162015E+00;
        lab_nodes_1D(2) =  - 0.867801053830347251000220202908E+00;
        lab_nodes_1D(3) =  - 0.728868599091326140584672400521E+00;
        lab_nodes_1D(4) =  - 0.550639402928647055316622705859E+00;
        lab_nodes_1D(5) =  - 0.342724013342712845043903403642E+00;
        lab_nodes_1D(6) =  - 0.116331868883703867658776709736E+00;
        lab_nodes_1D(7) =    0.116331868883703867658776709736E+00;
        lab_nodes_1D(8) =    0.342724013342712845043903403642E+00;
        lab_nodes_1D(9) =   0.550639402928647055316622705859E+00;
        lab_nodes_1D(10) =   0.728868599091326140584672400521E+00;
        lab_nodes_1D(11) =   0.867801053830347251000220202908E+00;
        lab_nodes_1D(12) =   0.959935045267260901355100162015E+00;
        lab_nodes_1D(13) =   1.0E+00;
    }

    else if (num == 15){
        lab_nodes_1D(0) =  - 1.0E+00;
        lab_nodes_1D(1) =  - 0.965245926503838572795851392070E+00;
        lab_nodes_1D(2) =  - 0.885082044222976298825401631482E+00;
        lab_nodes_1D(3) =  - 0.763519689951815200704118475976E+00;
        lab_nodes_1D(4) =  - 0.606253205469845711123529938637E+00;
        lab_nodes_1D(5) =  - 0.420638054713672480921896938739E+00;
        lab_nodes_1D(6) =  - 0.215353955363794238225679446273E+00;
        lab_nodes_1D(7) =    0.0E+00;
        lab_nodes_1D(8) =    0.215353955363794238225679446273E+00;
        lab_nodes_1D(9) =   0.420638054713672480921896938739E+00;
        lab_nodes_1D(10) =   0.606253205469845711123529938637E+00;
        lab_nodes_1D(11) =   0.763519689951815200704118475976E+00;
        lab_nodes_1D(12) =   0.885082044222976298825401631482E+00;
        lab_nodes_1D(13) =   0.965245926503838572795851392070E+00;
        lab_nodes_1D(14) =   1.0E+00;
    }

    else if (num == 16){
        lab_nodes_1D(0) =  - 1.0E+00;
        lab_nodes_1D(1) =  - 0.969568046270217932952242738367E+00;
        lab_nodes_1D(2) =  - 0.899200533093472092994628261520E+00;
        lab_nodes_1D(3) =  - 0.792008291861815063931088270963E+00;
        lab_nodes_1D(4) =  - 0.652388702882493089467883219641E+00;
        lab_nodes_1D(5) =  - 0.486059421887137611781890785847E+00;
        lab_nodes_1D(6) =  - 0.299830468900763208098353454722E+00;
        lab_nodes_1D(7) =  - 0.101326273521949447843033005046E+00;
        lab_nodes_1D(8) =    0.101326273521949447843033005046E+00;
        lab_nodes_1D(9) =   0.299830468900763208098353454722E+00;
        lab_nodes_1D(10) =   0.486059421887137611781890785847E+00;
        lab_nodes_1D(11) =   0.652388702882493089467883219641E+00;
        lab_nodes_1D(12) =   0.792008291861815063931088270963E+00;
        lab_nodes_1D(13) =   0.899200533093472092994628261520E+00;
        lab_nodes_1D(14) =   0.969568046270217932952242738367E+00;
        lab_nodes_1D(15) =   1.0E+00;
    }

    else if (num == 17){
        lab_nodes_1D(0) =  - 1.0E+00;
        lab_nodes_1D(1) =  - 0.973132176631418314156979501874E+00;
        lab_nodes_1D(2) =  - 0.910879995915573595623802506398E+00;
        lab_nodes_1D(3) =  - 0.815696251221770307106750553238E+00;
        lab_nodes_1D(4) =  - 0.691028980627684705394919357372E+00;
        lab_nodes_1D(5) =  - 0.541385399330101539123733407504E+00;
        lab_nodes_1D(6) =  - 0.372174433565477041907234680735E+00;
        lab_nodes_1D(7) =  - 0.189511973518317388304263014753E+00;
        lab_nodes_1D(8) =    0.0E+00;
        lab_nodes_1D(9) =   0.189511973518317388304263014753E+00;
        lab_nodes_1D(10) =   0.372174433565477041907234680735E+00;
        lab_nodes_1D(11) =   0.541385399330101539123733407504E+00;
        lab_nodes_1D(12) =   0.691028980627684705394919357372E+00;
        lab_nodes_1D(13) =   0.815696251221770307106750553238E+00;
        lab_nodes_1D(14) =   0.910879995915573595623802506398E+00;
        lab_nodes_1D(15) =   0.973132176631418314156979501874E+00;
        lab_nodes_1D(16) =   1.0E+00;
    }

    else if (num == 18){
        lab_nodes_1D(0) =  - 1.0E+00;
        lab_nodes_1D(1) =  - 0.976105557412198542864518924342E+00;
        lab_nodes_1D(2) =  - 0.920649185347533873837854625431E+00;
        lab_nodes_1D(3) =  - 0.835593535218090213713646362328E+00;
        lab_nodes_1D(4) =  - 0.723679329283242681306210365302E+00;
        lab_nodes_1D(5) =  - 0.588504834318661761173535893194E+00;
        lab_nodes_1D(6) =  - 0.434415036912123975342287136741E+00;
        lab_nodes_1D(7) =  - 0.266362652878280984167665332026E+00;
        lab_nodes_1D(8) =  - 0.897490934846521110226450100886E-01;
        lab_nodes_1D(9) =   0.897490934846521110226450100886E-01;
        lab_nodes_1D(10) =   0.266362652878280984167665332026E+00;
        lab_nodes_1D(11) =   0.434415036912123975342287136741E+00;
        lab_nodes_1D(12) =   0.588504834318661761173535893194E+00;
        lab_nodes_1D(13) =   0.723679329283242681306210365302E+00;
        lab_nodes_1D(14) =   0.835593535218090213713646362328E+00;
        lab_nodes_1D(15) =   0.920649185347533873837854625431E+00;
        lab_nodes_1D(16) =   0.976105557412198542864518924342E+00;
        lab_nodes_1D(17) =   1.0E+00;
    }

    else if (num == 19) {
        lab_nodes_1D(0)=   - 1.0E+00;
        lab_nodes_1D(1)=   - 0.978611766222080095152634063110E+00;
        lab_nodes_1D(2)=   - 0.928901528152586243717940258797E+00;
        lab_nodes_1D(3)=   - 0.852460577796646093085955970041E+00;
        lab_nodes_1D(4)=   - 0.751494202552613014163637489634E+00;
        lab_nodes_1D(5)=   - 0.628908137265220497766832306229E+00;
        lab_nodes_1D(6)=   - 0.488229285680713502777909637625E+00;
        lab_nodes_1D(7)=   - 0.333504847824498610298500103845E+00;
        lab_nodes_1D(8)=   - 0.169186023409281571375154153445E+00;
        lab_nodes_1D(9)=     0.0E+00;
        lab_nodes_1D(10) =   0.169186023409281571375154153445E+00;
        lab_nodes_1D(11) =   0.333504847824498610298500103845E+00;
        lab_nodes_1D(12) =   0.488229285680713502777909637625E+00;
        lab_nodes_1D(13) =   0.628908137265220497766832306229E+00;
        lab_nodes_1D(14) =   0.751494202552613014163637489634E+00;
        lab_nodes_1D(15) =   0.852460577796646093085955970041E+00;
        lab_nodes_1D(16) =   0.928901528152586243717940258797E+00;
        lab_nodes_1D(17) =   0.978611766222080095152634063110E+00;
        lab_nodes_1D(18) =   1.0E+00;

    } // end if

    else if (num == 20){
        lab_nodes_1D(0) =  - 1.0E+00;
        lab_nodes_1D(1) =  - 0.980743704893914171925446438584E+00;
        lab_nodes_1D(2) =  - 0.935934498812665435716181584931E+00;
        lab_nodes_1D(3) =  - 0.866877978089950141309847214616E+00;
        lab_nodes_1D(4) =  - 0.775368260952055870414317527595E+00;
        lab_nodes_1D(5) =  - 0.663776402290311289846403322971E+00;
        lab_nodes_1D(6) =  - 0.534992864031886261648135961829E+00;
        lab_nodes_1D(7) =  - 0.392353183713909299386474703816E+00;
        lab_nodes_1D(8) =  - 0.239551705922986495182401356927E+00;
        lab_nodes_1D(9) =  - 0.805459372388218379759445181596E-01;
        lab_nodes_1D(10) =   0.805459372388218379759445181596E-01;
        lab_nodes_1D(11) =   0.239551705922986495182401356927E+00;
        lab_nodes_1D(12) =   0.392353183713909299386474703816E+00;
        lab_nodes_1D(13) =   0.534992864031886261648135961829E+00;
        lab_nodes_1D(14) =   0.663776402290311289846403322971E+00;
        lab_nodes_1D(15) =   0.775368260952055870414317527595E+00;
        lab_nodes_1D(16) =   0.866877978089950141309847214616E+00;
        lab_nodes_1D(17) =   0.935934498812665435716181584931E+00;
        lab_nodes_1D(18) =   0.980743704893914171925446438584E+00;
        lab_nodes_1D(19) =   1.0E+00;
    }
} // end of labbato_nodes_1D function

// creates quadrature weights for Labatto polynomial
void labatto_weights_1D(
    c_array_t <real_t> &lab_weights_1D,  // Labbatto weights
    const int &num){                     // Interpolation order
    if (num == 1){
        lab_weights_1D(0) = 2.0;
    }
    else if (num == 2){
        lab_weights_1D(0) = 1.0;
        lab_weights_1D(1) = 1.0;
    }
    else if (num == 3){
        lab_weights_1D(0) = 1.0/3.0;
        lab_weights_1D(1) = 4.0/3.0;
        lab_weights_1D(2) = 1.0/3.0;
    }
    else if (num == 4){
        lab_weights_1D(0) = 1.0/6.0;
        lab_weights_1D(1) = 5.0/6.0;
        lab_weights_1D(2) = 5.0/6.0;
        lab_weights_1D(3) = 1.0/6.0;
    }
    else if (num == 5){
        lab_weights_1D(0) = 1.0/10.0;
        lab_weights_1D(1) = 49.0/90.0;
        lab_weights_1D(2) = 32.0/45.0;
        lab_weights_1D(3) = 49.0/90.0;
        lab_weights_1D(4) = 1.0/10.0;
    }
    else if (num == 6){
        lab_weights_1D(0) = 1.0/15.0;
        lab_weights_1D(1) = 1.0/30.0*(14.0 - sqrt(7.0));
        lab_weights_1D(2) = 1.0/30.0*(14.0 + sqrt(7.0));
        lab_weights_1D(3) = 1.0/30.0*(14.0 + sqrt(7.0));
        lab_weights_1D(4) = 1.0/30.0*(14.0 - sqrt(7.0));
        lab_weights_1D(5) = 1.0/15.0;
    }
    else if (num == 7){
        lab_weights_1D(0) =  0.476190476190476190476190476190E-01;
        lab_weights_1D(1) =  0.276826047361565948010700406290E+00;
        lab_weights_1D(2) =  0.431745381209862623417871022281E+00;
        lab_weights_1D(3) =  0.487619047619047619047619047619E+00;
        lab_weights_1D(4) =  0.431745381209862623417871022281E+00;
        lab_weights_1D(5) =  0.276826047361565948010700406290E+00;
        lab_weights_1D(6) =  0.476190476190476190476190476190E-01;
    }
    else if (num == 8){
        lab_weights_1D(0) =  0.357142857142857142857142857143E-01;
        lab_weights_1D(1) =  0.210704227143506039382991065776E+00;
        lab_weights_1D(2) =  0.341122692483504364764240677108E+00;
        lab_weights_1D(3) =  0.412458794658703881567052971402E+00;
        lab_weights_1D(4) =  0.412458794658703881567052971402E+00;
        lab_weights_1D(5) =  0.341122692483504364764240677108E+00;
        lab_weights_1D(6) =  0.210704227143506039382991065776E+00;
        lab_weights_1D(7) =  0.357142857142857142857142857143E-01;
    }
    else if (num == 9){
        lab_weights_1D(0) =  0.277777777777777777777777777778E-01;
        lab_weights_1D(1) =  0.165495361560805525046339720029E+00;
        lab_weights_1D(2) =  0.274538712500161735280705618579E+00;
        lab_weights_1D(3) =  0.346428510973046345115131532140E+00;
        lab_weights_1D(4) =  0.371519274376417233560090702948E+00;
        lab_weights_1D(5) =  0.346428510973046345115131532140E+00;
        lab_weights_1D(6) =  0.274538712500161735280705618579E+00;
        lab_weights_1D(7) =  0.165495361560805525046339720029E+00;
        lab_weights_1D(8) =  0.277777777777777777777777777778E-01;
    }
    else if (num == 10){
        lab_weights_1D(0) =  0.222222222222222222222222222222E-01;
        lab_weights_1D(1) =  0.133305990851070111126227170755E+00;
        lab_weights_1D(2) =  0.224889342063126452119457821731E+00;
        lab_weights_1D(3) =  0.292042683679683757875582257374E+00;
        lab_weights_1D(4) =  0.327539761183897456656510527917E+00;
        lab_weights_1D(5) =  0.327539761183897456656510527917E+00;
        lab_weights_1D(6) =  0.292042683679683757875582257374E+00;
        lab_weights_1D(7) =  0.224889342063126452119457821731E+00;
        lab_weights_1D(8) =  0.133305990851070111126227170755E+00;
        lab_weights_1D(9) =  0.222222222222222222222222222222E-01;
    }
    else if (num == 11){
        lab_weights_1D(0) =  0.181818181818181818181818181818E-01;
        lab_weights_1D(1) =  0.109612273266994864461403449580E+00;
        lab_weights_1D(2) =  0.187169881780305204108141521899E+00;
        lab_weights_1D(3) =  0.248048104264028314040084866422E+00;
        lab_weights_1D(4) =  0.286879124779008088679222403332E+00;
        lab_weights_1D(5) =  0.300217595455690693785931881170E+00;
        lab_weights_1D(6) =  0.286879124779008088679222403332E+00;
        lab_weights_1D(7) =  0.248048104264028314040084866422E+00;
        lab_weights_1D(8) =  0.187169881780305204108141521899E+00;
        lab_weights_1D(9) =  0.109612273266994864461403449580E+00;
        lab_weights_1D(10)=  0.181818181818181818181818181818E-01;
    }

    else if (num == 12){
        lab_weights_1D(0) =  0.151515151515151515151515151515E-01;
        lab_weights_1D(1) =  0.916845174131961306683425941341E-01;
        lab_weights_1D(2) =  0.157974705564370115164671062700E+00;
        lab_weights_1D(3) =  0.212508417761021145358302077367E+00;
        lab_weights_1D(4) =  0.251275603199201280293244412148E+00;
        lab_weights_1D(5) =  0.271405240910696177000288338500E+00;
        lab_weights_1D(6) =  0.271405240910696177000288338500E+00;
        lab_weights_1D(7) =  0.251275603199201280293244412148E+00;
        lab_weights_1D(8) =  0.212508417761021145358302077367E+00;
        lab_weights_1D(9) =  0.157974705564370115164671062700E+00;
        lab_weights_1D(10) = 0.916845174131961306683425941341E-01;
        lab_weights_1D(11) = 0.151515151515151515151515151515E-01;
    }

    else if (num == 13){
        lab_weights_1D(0) =  0.128205128205128205128205128205E-01;
        lab_weights_1D(1) =  0.778016867468189277935889883331E-01;
        lab_weights_1D(2) =  0.134981926689608349119914762589E+00;
        lab_weights_1D(3) =  0.183646865203550092007494258747E+00;
        lab_weights_1D(4) =  0.220767793566110086085534008379E+00;
        lab_weights_1D(5) =  0.244015790306676356458578148360E+00;
        lab_weights_1D(6) =  0.251930849333446736044138641541E+00;
        lab_weights_1D(7) =  0.244015790306676356458578148360E+00;
        lab_weights_1D(8) =  0.220767793566110086085534008379E+00;
        lab_weights_1D(9) =  0.183646865203550092007494258747E+00;
        lab_weights_1D(10) = 0.134981926689608349119914762589E+00;
        lab_weights_1D(11) = 0.778016867468189277935889883331E-01;
        lab_weights_1D(12) = 0.128205128205128205128205128205E-01;
    }

    else if (num == 14){
        lab_weights_1D(0) =  0.109890109890109890109890109890E-01;
        lab_weights_1D(1) =  0.668372844976812846340706607461E-01;
        lab_weights_1D(2) =  0.116586655898711651540996670655E+00;
        lab_weights_1D(3) =  0.160021851762952142412820997988E+00;
        lab_weights_1D(4) =  0.194826149373416118640331778376E+00;
        lab_weights_1D(5) =  0.219126253009770754871162523954E+00;
        lab_weights_1D(6) =  0.231612794468457058889628357293E+00;
        lab_weights_1D(7) =  0.231612794468457058889628357293E+00;
        lab_weights_1D(8) =  0.219126253009770754871162523954E+00;
        lab_weights_1D(9) =  0.194826149373416118640331778376E+00;
        lab_weights_1D(10) = 0.160021851762952142412820997988E+00;
        lab_weights_1D(11) = 0.116586655898711651540996670655E+00;
        lab_weights_1D(12) = 0.668372844976812846340706607461E-01;
        lab_weights_1D(13) = 0.109890109890109890109890109890E-01;
    }


    else if (num == 15){
        lab_weights_1D(0) =  0.952380952380952380952380952381E-02;
        lab_weights_1D(1) =  0.580298930286012490968805840253E-01;
        lab_weights_1D(2) =  0.101660070325718067603666170789E+00;
        lab_weights_1D(3) =  0.140511699802428109460446805644E+00;
        lab_weights_1D(4) =  0.172789647253600949052077099408E+00;
        lab_weights_1D(5) =  0.196987235964613356092500346507E+00;
        lab_weights_1D(6) =  0.211973585926820920127430076977E+00;
        lab_weights_1D(7) =  0.217048116348815649514950214251E+00;
        lab_weights_1D(8) =  0.211973585926820920127430076977E+00;
        lab_weights_1D(9) =  0.196987235964613356092500346507E+00;
        lab_weights_1D(10) = 0.172789647253600949052077099408E+00;
        lab_weights_1D(11) = 0.140511699802428109460446805644E+00;
        lab_weights_1D(12) = 0.101660070325718067603666170789E+00;
        lab_weights_1D(13) = 0.580298930286012490968805840253E-01;
        lab_weights_1D(14) = 0.952380952380952380952380952381E-02;
    }


    else if (num == 16){
        lab_weights_1D(0) =  0.833333333333333333333333333333E-02;
        lab_weights_1D(1) =  0.508503610059199054032449195655E-01;
        lab_weights_1D(2) =  0.893936973259308009910520801661E-01;
        lab_weights_1D(3) =  0.124255382132514098349536332657E+00;
        lab_weights_1D(4) =  0.154026980807164280815644940485E+00;
        lab_weights_1D(5) =  0.177491913391704125301075669528E+00;
        lab_weights_1D(6) =  0.193690023825203584316913598854E+00;
        lab_weights_1D(7) =  0.201958308178229871489199125411E+00;
        lab_weights_1D(8) =  0.201958308178229871489199125411E+00;
        lab_weights_1D(9) =  0.193690023825203584316913598854E+00;
        lab_weights_1D(10) = 0.177491913391704125301075669528E+00;
        lab_weights_1D(11) = 0.154026980807164280815644940485E+00;
        lab_weights_1D(12) = 0.124255382132514098349536332657E+00;
        lab_weights_1D(13) = 0.893936973259308009910520801661E-01;
        lab_weights_1D(14) = 0.508503610059199054032449195655E-01;
        lab_weights_1D(15) = 0.833333333333333333333333333333E-02;
    }


    else if (num == 17){
        lab_weights_1D(0) =  0.735294117647058823529411764706E-02;
        lab_weights_1D(1) =  0.449219405432542096474009546232E-01;
        lab_weights_1D(2) =  0.791982705036871191902644299528E-01;
        lab_weights_1D(3) =  0.110592909007028161375772705220E+00;
        lab_weights_1D(4) =  0.137987746201926559056201574954E+00;
        lab_weights_1D(5) =  0.160394661997621539516328365865E+00;
        lab_weights_1D(6) =  0.177004253515657870436945745363E+00;
        lab_weights_1D(7) =  0.187216339677619235892088482861E+00;
        lab_weights_1D(8) =  0.190661874753469433299407247028E+00;
        lab_weights_1D(9) =  0.187216339677619235892088482861E+00;
        lab_weights_1D(10) = 0.177004253515657870436945745363E+00;
        lab_weights_1D(11) = 0.160394661997621539516328365865E+00;
        lab_weights_1D(12) = 0.137987746201926559056201574954E+00;
        lab_weights_1D(13) = 0.110592909007028161375772705220E+00;
        lab_weights_1D(14) = 0.791982705036871191902644299528E-01;
        lab_weights_1D(15) = 0.449219405432542096474009546232E-01;
        lab_weights_1D(16) = 0.735294117647058823529411764706E-02;
    }

    else if (num == 18){
        lab_weights_1D(0) =  0.653594771241830065359477124183E-02;
        lab_weights_1D(1) =  0.399706288109140661375991764101E-01;
        lab_weights_1D(2) =  0.706371668856336649992229601678E-01;
        lab_weights_1D(3) =  0.990162717175028023944236053187E-01;
        lab_weights_1D(4) =  0.124210533132967100263396358897E+00;
        lab_weights_1D(5) =  0.145411961573802267983003210494E+00;
        lab_weights_1D(6) =  0.161939517237602489264326706700E+00;
        lab_weights_1D(7) =  0.173262109489456226010614403827E+00;
        lab_weights_1D(8) =  0.179015863439703082293818806944E+00;
        lab_weights_1D(9) =  0.179015863439703082293818806944E+00;
        lab_weights_1D(10) = 0.173262109489456226010614403827E+00;
        lab_weights_1D(11) = 0.161939517237602489264326706700E+00;
        lab_weights_1D(12) = 0.145411961573802267983003210494E+00;
        lab_weights_1D(13) = 0.124210533132967100263396358897E+00;
        lab_weights_1D(14) = 0.990162717175028023944236053187E-01;
        lab_weights_1D(15) = 0.706371668856336649992229601678E-01;
        lab_weights_1D(16) = 0.399706288109140661375991764101E-01;
        lab_weights_1D(17) = 0.653594771241830065359477124183E-02;
    }

    else if (num == 19) {
        lab_weights_1D(0) =  0.584795321637426900584795321637E-02;
        lab_weights_1D(1) =  0.357933651861764771154255690351E-01;
        lab_weights_1D(2) =  0.633818917626297368516956904183E-01;
        lab_weights_1D(3) =  0.891317570992070844480087905562E-01;
        lab_weights_1D(4) =  0.112315341477305044070910015464E+00;
        lab_weights_1D(5) =  0.132267280448750776926046733910E+00;
        lab_weights_1D(6) =  0.148413942595938885009680643668E+00;
        lab_weights_1D(7) =  0.160290924044061241979910968184E+00;
        lab_weights_1D(8) =  0.167556584527142867270137277740E+00;
        lab_weights_1D(9) =  0.170001919284827234644672715617E+00;
        lab_weights_1D(10) = 0.167556584527142867270137277740E+00;
        lab_weights_1D(11) = 0.160290924044061241979910968184E+00;
        lab_weights_1D(12) = 0.148413942595938885009680643668E+00;
        lab_weights_1D(13) = 0.132267280448750776926046733910E+00;
        lab_weights_1D(14) = 0.112315341477305044070910015464E+00;
        lab_weights_1D(15) = 0.891317570992070844480087905562E-01;
        lab_weights_1D(16) = 0.633818917626297368516956904183E-01;
        lab_weights_1D(17) = 0.357933651861764771154255690351E-01;
        lab_weights_1D(18) = 0.584795321637426900584795321637E-02;
    } // end if

    else if (num == 20) {
        lab_weights_1D(0) =  0.526315789473684210526315789474E-02;
        lab_weights_1D(1) =  0.322371231884889414916050281173E-01;
        lab_weights_1D(2) =  0.571818021275668260047536271732E-01;
        lab_weights_1D(3) =  0.806317639961196031447768461137E-01;
        lab_weights_1D(4) =  0.101991499699450815683781205733E+00;
        lab_weights_1D(5) =  0.120709227628674725099429705002E+00;
        lab_weights_1D(6) =  0.136300482358724184489780792989E+00;
        lab_weights_1D(7) =  0.148361554070916825814713013734E+00;
        lab_weights_1D(8) =  0.156580102647475487158169896794E+00;
        lab_weights_1D(9) =  0.160743286387845749007726726449E+00;
        lab_weights_1D(10) = 0.160743286387845749007726726449E+00;
        lab_weights_1D(11) = 0.156580102647475487158169896794E+00;
        lab_weights_1D(12) = 0.148361554070916825814713013734E+00;
        lab_weights_1D(13) = 0.136300482358724184489780792989E+00;
        lab_weights_1D(14) = 0.120709227628674725099429705002E+00;
        lab_weights_1D(15) = 0.101991499699450815683781205733E+00;
        lab_weights_1D(16) = 0.806317639961196031447768461137E-01;
        lab_weights_1D(17) = 0.571818021275668260047536271732E-01;
        lab_weights_1D(18) = 0.322371231884889414916050281173E-01;
        lab_weights_1D(19) = 0.526315789473684210526315789474E-02;
    } // end if
} // end of labatto_weights_1D function

// Create set of quadrature weights for a line through the subcells
void length_weights(
    c_array_t <real_t> &len_weights_1D,  // Labbatto weights
    c_array_t <real_t> &lab_weights_1D,  // Labbatto weights
    c_array_t <real_t> &lab_nodes_1D,
    const int &p_order){

    real_t alpha1 = (lab_nodes_1D(1) - lab_nodes_1D(0) - lab_weights_1D(0))
                    /lab_weights_1D(1);

    real_t alpha2 = (lab_nodes_1D(2) - lab_nodes_1D(0) 
                  - (lab_weights_1D(0)+ lab_weights_1D(1)))
                  /  lab_weights_1D(2);


    if(p_order == 1){
        len_weights_1D(0) = 1.0;
        len_weights_1D(1) = 1.0;
    }


    if(p_order == 2){
        len_weights_1D(0) = lab_weights_1D(0) + alpha1*lab_weights_1D(1);
        len_weights_1D(1) = (1.0-alpha1)*lab_weights_1D(1) + 0.5*lab_weights_1D(2);
        len_weights_1D(2) = len_weights_1D(1);
        len_weights_1D(3) = len_weights_1D(0);
    }

    if(p_order == 3){
        len_weights_1D(0) = lab_weights_1D(0) + alpha1*lab_weights_1D(1);
        len_weights_1D(1) = (1.0 - alpha1)*lab_weights_1D(1) + alpha2*(lab_weights_1D(2));
        len_weights_1D(2) = (1.0 - alpha2)*lab_weights_1D(2) + 0.5*lab_weights_1D(3);
        len_weights_1D(3) = len_weights_1D(2);
        len_weights_1D(4) = len_weights_1D(1);
        len_weights_1D(5) = len_weights_1D(0);
    }
}

// Create set of quadrature weights for the subcells
void sub_weights(
    c_array_t <real_t> &sub_weights_1D,  // Labbatto weights
    c_array_t <real_t> &lab_weights_1D,  // Labbatto weights
    c_array_t <real_t> &lab_nodes_1D,
    const int &p_order){

    real_t alpha1 = (lab_nodes_1D(1) - lab_nodes_1D(0) - lab_weights_1D(0))
                    /lab_weights_1D(1);

    real_t alpha2 = (lab_nodes_1D(2) - lab_nodes_1D(0) 
                  - (lab_weights_1D(0)+ lab_weights_1D(1)))
                  /  lab_weights_1D(2);


    if(p_order == 1){
        sub_weights_1D(0) = lab_weights_1D(0);
        sub_weights_1D(1) = 0.5*lab_weights_1D(1);
        sub_weights_1D(2) = sub_weights_1D(1);
        sub_weights_1D(3) = sub_weights_1D(0);
    }


    if(p_order == 2){
        sub_weights_1D(0) = lab_weights_1D(0);
        sub_weights_1D(1) = alpha1*lab_weights_1D(1);
        sub_weights_1D(2) = (1.0-alpha1)*lab_weights_1D(1);
        sub_weights_1D(3) = 0.5*lab_weights_1D(2);
        sub_weights_1D(4) = sub_weights_1D(3);
        sub_weights_1D(5) = sub_weights_1D(2);
        sub_weights_1D(6) = sub_weights_1D(1);
        sub_weights_1D(7) = sub_weights_1D(0);
    }

    if(p_order == 3){
        sub_weights_1D(0) = lab_weights_1D(0);
        sub_weights_1D(1) = alpha1*lab_weights_1D(1);
        sub_weights_1D(2) = (1.0-alpha1)*lab_weights_1D(1);
        sub_weights_1D(3) = alpha2*lab_weights_1D(2);
        sub_weights_1D(4) = (1.0-alpha2)*lab_weights_1D(2);
        sub_weights_1D(5) = 0.5*lab_weights_1D(3);
        sub_weights_1D(6) = sub_weights_1D(5);
        sub_weights_1D(7) = sub_weights_1D(4);
        sub_weights_1D(8) = sub_weights_1D(3);
        sub_weights_1D(9) = sub_weights_1D(2);
        sub_weights_1D(10) = sub_weights_1D(1);
        sub_weights_1D(11) = sub_weights_1D(0);
    }
}

// Some linear algebra snippets
void mat_inverse(
    c_array_t <real_t> &mat_inv,
    c_array_t <real_t> &matrix){

    double A_11 = matrix(1, 1)*matrix(2, 2) 
                - matrix(2, 1)*matrix(1, 2);

    double A_22 = matrix(2, 2)*matrix(0, 0) 
                - matrix(0, 2)*matrix(2, 0);

    double A_33 = matrix(0, 0)*matrix(1, 1) 
                - matrix(1, 0)*matrix(0, 1);

    double A_12 = matrix(2, 1)*matrix(0, 2) 
                - matrix(0, 1)*matrix(2, 2);

    double A_23 = matrix(0, 2)*matrix(1, 0) 
                - matrix(1, 2)*matrix(0, 0);

    double A_31 = matrix(1, 0)*matrix(2, 1) 
                - matrix(2, 0)*matrix(1, 1);

    double A_21 = matrix(1, 2)*matrix(2, 0) 
                - matrix(1, 0)*matrix(2, 2);

    double A_32 = matrix(2, 0)*matrix(0, 1) 
                - matrix(2, 1)*matrix(0, 0);

    double A_13 = matrix(0, 1)*matrix(1, 1) 
                - matrix(0, 2)*matrix(1, 1);

    double  det = matrix(0, 0)*A_11 + matrix(1, 0)*A_21 
                + matrix(2, 0)*A_31;

    mat_inv(0, 0) = A_11/det;
    mat_inv(0, 1) = A_12/det;
    mat_inv(0, 2) = A_13/det; 
    mat_inv(1, 0) = A_21/det;
    mat_inv(1, 1) = A_22/det;
    mat_inv(1, 2) = A_23/det;
    mat_inv(2, 0) = A_31/det;
    mat_inv(2, 1) = A_32/det;
    mat_inv(2, 2) = A_33/det;
}

void mat_mult(
    c_array_t <real_t> &result,
    c_array_t <real_t> &A,
    c_array_t <real_t> &B){

    result(0, 0) = A(0, 0)*B(0, 0) + A(0, 1)*B(1, 0) + A(0, 2)*B(2, 0);
    result(0, 1) = A(0, 0)*B(0, 1) + A(0, 1)*B(1, 1) + A(0, 2)*B(2, 1);
    result(0, 2) = A(0, 0)*B(0, 2) + A(0, 1)*B(1, 2) + A(0, 2)*B(2, 2);
    result(1, 0) = A(1, 0)*B(0, 0) + A(1, 1)*B(1, 0) + A(1, 2)*B(2, 0);
    result(1, 1) = A(1, 0)*B(0, 1) + A(1, 1)*B(1, 1) + A(1, 2)*B(2, 1);
    result(1, 2) = A(1, 0)*B(0, 2) + A(1, 1)*B(1, 2) + A(1, 2)*B(2, 2);
    result(2, 0) = A(2, 0)*B(0, 0) + A(2, 1)*B(1, 0) + A(2, 2)*B(2, 0);
    result(2, 1) = A(2, 0)*B(0, 1) + A(2, 1)*B(1, 1) + A(2, 2)*B(2, 1);
    result(2, 2) = A(2, 0)*B(0, 2) + A(2, 1)*B(1, 2) + A(2, 2)*B(2, 2);
}

void mat_trans(
    c_array_t <real_t> &trans,
    c_array_t <real_t> &mat){

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            trans(i, j) = mat(j, i);    
        }
    }
}

void set_nodes_wgts(
    c_array_t <real_t> &lab_nodes_1D,
    c_array_t <real_t> &lab_weights_1D,
    c_array_t <real_t> &len_weights_1D,
    c_array_t <real_t> &sub_weights_1D, 
    const int p_order){

    int num_g_pts_1d = 2 * p_order + 1;

    // Gauss Lobatto  Weights
    labatto_weights_1D(lab_weights_1D, num_g_pts_1d);

    // Gauss Lobatto  Nodes
    labatto_nodes_1D(lab_nodes_1D, num_g_pts_1d);

    // Get the lengths of each edge
    length_weights(len_weights_1D, lab_weights_1D, lab_nodes_1D, p_order);

    // Get the lengths of the partitioned edges (sub_weights)
    sub_weights(sub_weights_1D, lab_weights_1D, lab_nodes_1D, p_order);
}

void refine_mesh(
    mesh_t& init_mesh, 
    mesh_t& mesh, 
    const int p_order, 
    const int rk_num_bins,
    const int dim){

    // High order mesh parameters
    int num_sub_1d;         // num subcells in 1d
    int num_g_pts_1d;       // num gauss points in 1d

    int num_g_pts;            // number of gauss points
    int num_subcells_per_elem; // number subcells in an element

    int num_elem = init_mesh.num_cells();

    if(p_order == 0){
                
        num_sub_1d   = 1;
        num_g_pts_1d = 2;
        
        num_g_pts  = pow(num_g_pts_1d, dim);
        num_subcells_per_elem = pow((num_g_pts_1d-1), dim);
    }

    else{

        num_sub_1d = p_order*2;         // num subcells in 1d
        num_g_pts_1d = 2 * p_order + 1; // num gauss points in 1d

        num_g_pts  = pow(num_g_pts_1d, dim);            // number of gauss points
        num_subcells_per_elem = pow((num_sub_1d), dim); // number subcells in an element
    }

    //  ---------------------------------------------------------------------------
    //  Initailize Element and cell information in on high order mesh
    //  ---------------------------------------------------------------------------

    mesh.init_element(p_order, dim, num_elem, rk_num_bins);
    mesh.init_cells(num_elem*num_subcells_per_elem, rk_num_bins);
    mesh.init_gauss_pts();

    //  ---------------------------------------------------------------------------
    //  Generate point positiont in reference space to map onto initial mesh
    //  ---------------------------------------------------------------------------

    auto temp_pts = c_array_t<real_t> (num_g_pts_1d, num_g_pts_1d, num_g_pts_1d, 3);

    double dx = 2.0/((double)num_g_pts_1d - 1.0);  // len/(num_nodes-1)
    double dy = 2.0/((double)num_g_pts_1d - 1.0);  // len/(num_nodes-1)
    double dz = 2.0/((double)num_g_pts_1d - 1.0);  // len/(num_nodes-1)

    for(int k = 0; k < num_g_pts_1d; k++){
        for(int j = 0; j < num_g_pts_1d; j++){
            for(int i = 0; i < num_g_pts_1d; i++){
                temp_pts(i, j, k, 0) = -1.0 + (double)i*dx;
                temp_pts(i, j, k, 1) = -1.0 + (double)j*dy;
                temp_pts(i, j, k, 2) = -1.0 + (double)k*dz;
            }
        }
    }


    //  ---------------------------------------------------------------------------
    //  Map new points to real space using basis functions
    //  ---------------------------------------------------------------------------

    // temp array to hold positions
    real_t * temp_gauss_point_coords; 
    temp_gauss_point_coords = new real_t[num_elem*num_g_pts*dim];
    auto g_points_in_mesh = view_c_array <real_t> (temp_gauss_point_coords, num_elem*num_g_pts, dim);


    // Reference node positions for element (currently p1, replace with element library)
    real_t ref_vert[8][3] = // listed as {Xi, Eta, Mu}
        {
        // Bottom Nodes
        {-1.0, -1.0, -1.0},// 0
        {+1.0, -1.0, -1.0},// 1
        {-1.0, +1.0, -1.0},// 2
        {+1.0, +1.0, -1.0},// 3
        // Top Nodes
        {-1.0, -1.0, +1.0},// 4
        {+1.0, -1.0, +1.0},// 5
        {-1.0, +1.0, +1.0},// 6
        {+1.0, +1.0, +1.0},// 7
        };


    real_t basis[8];    // basis function evaluation value


    // Inital mesh node positions
    real_t x_init_a[8];
    real_t y_init_a[8];
    real_t z_init_a[8];

    auto x_init = view_c_array <real_t> (x_init_a, 8);
    auto y_init = view_c_array <real_t> (y_init_a, 8);
    auto z_init = view_c_array <real_t> (z_init_a, 8);

    int g_point_count = 0; 

    // mapping points using the basis functions
    for(int elem_gid = 0; elem_gid < num_elem; elem_gid++){

        // Assign initial positions from initial mesh
        for(int node_lid = 0; node_lid < 8; node_lid++){
            x_init(node_lid) = init_mesh.node_coords(0, init_mesh.nodes_in_cell(elem_gid, node_lid), 0);
            y_init(node_lid) = init_mesh.node_coords(0, init_mesh.nodes_in_cell(elem_gid, node_lid), 1);
            z_init(node_lid) = init_mesh.node_coords(0, init_mesh.nodes_in_cell(elem_gid, node_lid), 2);
        }

        // Walk over gauss points as i,j,k mesh to calculate basis 
        for(int k = 0; k < num_g_pts_1d; k++){
            for(int j = 0; j < num_g_pts_1d; j++){
                for(int i = 0; i < num_g_pts_1d; i++){
                    
                    for (int vert_lid = 0; vert_lid < 8; vert_lid++ ){
                        basis[vert_lid] = 1.0/8.0
                                        * (1.0 + temp_pts(i, j, k, 0)*ref_vert[vert_lid][0])
                                        * (1.0 + temp_pts(i, j, k, 1)*ref_vert[vert_lid][1])
                                        * (1.0 + temp_pts(i, j, k, 2)*ref_vert[vert_lid][2]);
                    }

                    g_points_in_mesh(g_point_count, 0) = 0.0;
                    g_points_in_mesh(g_point_count, 1) = 0.0;
                    g_points_in_mesh(g_point_count, 2) = 0.0;

                    // Walk over vertices to map new points onto mesh
                    for (int vert_lid = 0; vert_lid < 8; vert_lid++ ){
                        g_points_in_mesh(g_point_count, 0) += basis[vert_lid]*x_init(vert_lid);
                        g_points_in_mesh(g_point_count, 1) += basis[vert_lid]*y_init(vert_lid);
                        g_points_in_mesh(g_point_count, 2) += basis[vert_lid]*z_init(vert_lid);
                    } // end for vert_lid

                    g_point_count++;  
                }
            }
        }
    }

    //  ---------------------------------------------------------------------------
    //  Hash x, y, and x coordinates to eliminate double counted points for 
    //  node index. 
    //  ---------------------------------------------------------------------------


    real_t pos_max[dim];
    real_t pos_min[dim];

    for(int i = 0; i < dim; i++){
        pos_max[i] = -1.0E16;
        pos_min[i] =  1.0E16;
    }

    real_t position[3]; 

    // Get min and max points in the mesh
    for(int point = 0; point< init_mesh.num_nodes(); point++){
        for(int i = 0; i < dim; i++){
            
            position[i] = init_mesh.node_coords(0, point, i);
            pos_max[i] = fmax(pos_max[i], position[i]);
            pos_min[i] = fmin(pos_min[i], position[i]);
        }
    }

    // get minimum distance between any two points (WARNING: ONLY WORKS IN 3D)
    real_t dist_min;
    real_t dist_max;
    real_t cell_nodes[24];
    
    auto vert1 = view_c_array <real_t> (cell_nodes, 8, 3);
    
    real_t distance[28]; 
    auto dist = view_c_array <real_t> (distance, 28);

    for (int cell_gid = 0; cell_gid < init_mesh.num_cells(); cell_gid++){
        
        // Getting the coordinates of the element
        for(int node = 0; node < 8; node++){
            for (int dim = 0; dim < 3; dim++)
                vert1(node, dim) = init_mesh.node_coords(0, init_mesh.nodes_in_cell(cell_gid, node), dim);
        }

        // loop conditions needed for distance calculation
        int countA = 0;
        int countB = 1;
        int a;
        int b;
        int loop = 0;
        
        
        // Solving for the magnitude of distance between each node
        for (int i = 0; i < 28; i++){
            
            a = countA;
            b = countB;
            
            // returns magnitude of distance between each node, 28 total options
            dist(i) = fabs(sqrt(( pow((vert1(b,0) - vert1(a,0)), 2.0)
                                + pow((vert1(b,1) - vert1(a,1)), 2.0)
                                + pow((vert1(b,2) - vert1(a,2)), 2.0))));

            countB++;
            countA++;
            
            //tricky indexing
            if (countB > 7) {
                loop++;
                countB = 1 + loop;
                countA = 0;
            }
        }

        dist_min = 1e64;
        dist_max = 0.0;
        
        for(int i = 0; i < 28; ++i){
            dist_min = fmin(dist(i),dist_min);
            dist_max = fmax(dist(i),dist_max);
        }
    }

    std::cout<<"Min_dist = "<<dist_min<<std::endl;

    // Number of subcells per dimension to be created
    real_t sub;
    
    if (p_order == 0) sub = 1.0;
    else sub = p_order;

    real_t h = dist_min/(3.0*(sub)); // number of subdivisions between any two points

    // Define number of bins in each direction
    real_t float_bins[3];
    for(int i = 0; i < dim; i++){

        float_bins[i] = fmax(1e-16, (pos_max[i] - pos_min[i] + 1.0 + 1e-14)/h); //1e-14
    }

    // Convert # of bins to ints
    int int_bins[3];
    for(int i = 0; i < dim; i++){

        int_bins[i] = (int)float_bins[i];
    }

    real_t float_idx[3];   // float values for index
    int int_idx[3];        // int values for index
    
    int key;            // hash key 
    int max_key = 0;    // larges hash key value


    // Getting hash keys from x,y,z positions and largest key value
    int h_keys[num_g_pts*num_elem];
    for(int g_pt = 0; g_pt < num_g_pts*num_elem; g_pt++){
        
        real_t coords[3];
        for(int i = 0; i < dim; i++){

            coords[i]    = g_points_in_mesh(g_pt, i);
            float_idx[i] = fmax(1e-16, (coords[i] - pos_min[i] + 1e-14)/(h));
            int_idx[i] = (int)float_idx[i];
        }

        
        // i + j*num_x + k*num_x*num_y
        if (dim == 2){
            key = int_idx[0] + int_idx[1]*int_bins[0];
        }

        else{
            key = int_idx[0] 
                + int_idx[1]*int_bins[0] 
                + int_idx[2]*int_bins[0]*int_bins[1];
        }

        h_keys[g_pt] = key;
        max_key = std::max(max_key, key);
    }

    // Allocating array for hash table
    int * hash; 
    hash = new int[max_key+10]; 

    // Initializing values at key positions to zero
    for(int g_pt = 0; g_pt < num_g_pts*num_elem; g_pt++){
        hash[h_keys[g_pt]] = 0;
    }

    // Temporary array for gauss->node map and node->gauss map

    int * node_to_gauss_map;

    node_to_gauss_map = new int[num_g_pts*num_elem];
    // gauss_node_map    = new int[num_g_pts*num_elem];

    // counters
    int num_nodes = 0;
    int node_gid = 0;

    // walk over all gauss points 
    for(int g_pt = 0; g_pt < num_g_pts*num_elem; g_pt++){
        
        // Subtract 1 every time the index is touched
        if(hash[h_keys[g_pt]] <= 0){
            hash[h_keys[g_pt]] += -1;
        }
        
        // If this is the first time the index is touched add to 
        // node_to_gauss_map (WARNING: ONLY THE FIRST TIME IS COUNTED)
        // and index the number of nodes 
        if(hash[h_keys[g_pt]] == -1){

            node_to_gauss_map[num_nodes] = g_pt;
            num_nodes++;
        }

        // If this index has been touched before, replace hash value with
        // node id
        if(hash[h_keys[g_pt]] <= -1){
            
            hash[h_keys[g_pt]] = node_gid;
            
            // gauss_node_map[g_pt] = node_gid;
            mesh.node_in_gauss(g_pt) = node_gid;

            node_gid++;
        }

        // If hash value is positive, then the value is the index
        // for the single node assiciated with this g_point
        else{
            // gauss_node_map[g_pt] = hash[h_keys[g_pt]];
            mesh.node_in_gauss(g_pt) = hash[h_keys[g_pt]];
        }
    }

    // remove hash table
    delete[] hash;

    // Initialize nodes on sub_mesh
    mesh.init_nodes(num_nodes, rk_num_bins);

    //  ---------------------------------------------------------------------------
    //  Write position to nodes 
    //  ---------------------------------------------------------------------------
    
    for(int node_gid=0; node_gid<num_nodes; node_gid++){
        for(int i = 0; i < dim; i++){

            mesh.node_coords(0, node_gid, i) = g_points_in_mesh(node_to_gauss_map[node_gid], i);
        }
    }

    // delete unneeded arrays
    delete[] temp_gauss_point_coords;
    delete[] node_to_gauss_map;



    //  ---------------------------------------------------------------------------
    //  Get gauss points and nodes associated with each cell, 
    //  as well as the cells associated with each element
    //  ---------------------------------------------------------------------------

    // auto gauss_id_in_cell = c_array_t<int> (sub_mesh.num_cells(), num_sub_1d*num_sub_1d*num_sub_1d, 8);
    int sub_in_elem = num_sub_1d*num_sub_1d*num_sub_1d;
    int * gauss_id_in_cell;

    gauss_id_in_cell = new int[mesh.num_cells()*sub_in_elem*8];
    auto gauss_in_cell = view_c_array<int> (gauss_id_in_cell, mesh.num_cells(), sub_in_elem, 8);

    int p0, p1, p2, p3, p4, p5, p6, p7;
    p0 = p1 = p2 = p3 = p4 = p5 = p6 = p7 = 0;
    
    int num_1d = num_g_pts_1d;
    int cell_index = 0;
    int cell_mesh_index = 0;

    for(int elem_gid = 0; elem_gid < num_elem; elem_gid++){
        for(int k = 0; k < num_sub_1d; k++){
            for(int j = 0; j < num_sub_1d; j++){
                for(int i = 0; i < num_sub_1d; i++){

                    // The p# point to a global gauss point index before double counting 
                    p0 = (i)     + (j)*num_1d   + (k)*num_1d*num_1d;
                    p1 = (i+1)   + (j)*num_1d   + (k)*num_1d*num_1d;
                    p2 = (i)     + (j+1)*num_1d + (k)*num_1d*num_1d;
                    p3 = (i+1)   + (j+1)*num_1d + (k)*num_1d*num_1d;
                    p4 = (i)     + (j)*num_1d   + (k+1)*num_1d*num_1d;
                    p5 = (i+1)   + (j)*num_1d   + (k+1)*num_1d*num_1d;
                    p6 = (i)     + (j+1)*num_1d + (k+1)*num_1d*num_1d;
                    p7 = (i+1)   + (j+1)*num_1d + (k+1)*num_1d*num_1d;

                    p0 += num_1d*num_1d*num_1d*(elem_gid); 
                    p1 += num_1d*num_1d*num_1d*(elem_gid); 
                    p2 += num_1d*num_1d*num_1d*(elem_gid); 
                    p3 += num_1d*num_1d*num_1d*(elem_gid); 
                    p4 += num_1d*num_1d*num_1d*(elem_gid); 
                    p5 += num_1d*num_1d*num_1d*(elem_gid); 
                    p6 += num_1d*num_1d*num_1d*(elem_gid); 
                    p7 += num_1d*num_1d*num_1d*(elem_gid); 

                    cell_index = i + j*num_sub_1d + k*num_sub_1d*num_sub_1d;

                    cell_mesh_index = cell_index + num_sub_1d*num_sub_1d*num_sub_1d*(elem_gid);

                    // if(cell_mesh_index != elem_gid) std::cout<<"ERROR IN REFINE MESH"<<std::endl;

                    gauss_in_cell(elem_gid, cell_index, 0) = p0;
                    gauss_in_cell(elem_gid, cell_index, 1) = p1;
                    gauss_in_cell(elem_gid, cell_index, 2) = p2;
                    gauss_in_cell(elem_gid, cell_index, 3) = p3;
                    gauss_in_cell(elem_gid, cell_index, 4) = p4;
                    gauss_in_cell(elem_gid, cell_index, 5) = p5;
                    gauss_in_cell(elem_gid, cell_index, 6) = p6;
                    gauss_in_cell(elem_gid, cell_index, 7) = p7;

                    mesh.cells_in_elem(elem_gid, cell_index) = cell_mesh_index;

                    mesh.elems_in_cell(cell_mesh_index) = elem_gid;

                }
            }
        }
    }



    int cell_gid = 0;
    int p[8];
    // for each cell read the list of associated nodes
    for(int elem = 0; elem < num_elem; elem++){
        for(int k = 0; k < num_sub_1d; k++){
            for(int j = 0; j < num_sub_1d; j++){
                for(int i = 0; i < num_sub_1d; i++){


                    // The p# point to a global gauss point index before double counting 
                    p[0] = (i)     + (j)*num_1d   + (k)*num_1d*num_1d;
                    p[1] = (i+1)   + (j)*num_1d   + (k)*num_1d*num_1d;
                    p[2] = (i)     + (j+1)*num_1d + (k)*num_1d*num_1d;
                    p[3] = (i+1)   + (j+1)*num_1d + (k)*num_1d*num_1d;
                    p[4] = (i)     + (j)*num_1d   + (k+1)*num_1d*num_1d;
                    p[5] = (i+1)   + (j)*num_1d   + (k+1)*num_1d*num_1d;
                    p[6] = (i)     + (j+1)*num_1d + (k+1)*num_1d*num_1d;
                    p[7] = (i+1)   + (j+1)*num_1d + (k+1)*num_1d*num_1d;


                    for (int idx = 0; idx < 8; idx++){
                        p[idx] += num_1d*num_1d*num_1d*(elem); 
                    }

                    for (int node_lid = 0; node_lid < 8; node_lid++){
                        mesh.nodes_in_cell(cell_gid, node_lid) = mesh.node_in_gauss(p[node_lid]);
                    }
                    // incriment global index for cell
                    cell_gid++;
                }
            }
        }
    }

    delete [] gauss_id_in_cell;

    mesh.build_connectity();
}


namespace swage{

// Used by Gauss 2/3/4D to set quadrature points
void line_gauss_info(
    real_t &x, 
    real_t &w, 
    int &m, 
    int &p){

    m--;
    // creating matrices for weights and points
    real_t g2[2],w3[3],g3[3],w4[4],g4[4],g5[5],w5[5], 
      g6[6],w6[6],g7[7],w7[7],g8[8],w8[8];

    g2[0] = -1.0/sqrt(3.0);
    g2[1] = -1.0*g2[0]; 

    w3[0] = 5.0/9.0;
    w3[1] = 8.0/9.0;
    w3[2] = w3[0]; 

    g3[0] = -sqrt(3.0/5.0);
    g3[1] = 0.0;
    g3[2] = -1.0*g3[0]; 

    w4[0] = (1./2.)-sqrt(5./6.)/6.;
    w4[1] = (1./2.)+sqrt(5./6.)/6.;
    w4[2] = w4[1]; 
    w4[3] = w4[0]; 

    g4[0] = -sqrt((3.+2.*sqrt(6./5.))/7.);
    g4[1] = -sqrt((3.-2.*sqrt(6./5.))/7.);
    g4[2] = -g4[1]; 
    g4[3] = -g4[0]; 

    w5[0] = (322.-13.*sqrt(70.))/900.;
    w5[1] = (322.+13.*sqrt(70.))/900.;
    w5[2] =  512./900.;
    w5[3] = w5[1]; 
    w5[4] = w5[0]; 

    g5[0] = -sqrt(5.+2.*sqrt(10./7.))/3.;
    g5[1] = -sqrt(5.-2.*sqrt(10./7.))/3.;
    g5[2] = 0.0;
    g5[3] = -g5[1]; 
    g5[4] = -g5[0]; 

    w6[0] = 0.1713244923791704;
    w6[1] = 0.3607615730481386;
    w6[2] = 0.4679139345726910; 
    w6[3] = w6[2];
    w6[4] = w6[1];
    w6[5] = w6[0];

    g6[0] = -0.9324695142031521;
    g6[1] = -0.6612093864662645;
    g6[2] = -0.2386191860831969;
    g6[3] = -g6[2];
    g6[4] = -g6[1];
    g6[5] = -g6[0];

    w7[0] = 0.1294849661688697;
    w7[1] = 0.2797053914892766;
    w7[2] = 0.3818300505051189;
    w7[3] = 0.4179591836734694;
    w7[4] = w7[2];
    w7[5] = w7[1];
    w7[6] = w7[0];

    g7[0] = -0.9491079123427585;
    g7[1] = -0.7415311855993945;
    g7[2] = -0.4058451513773972;
    g7[3] =  0.0;
    g7[4] = -g7[2];
    g7[5] = -g7[1];
    g7[6] = -g7[0];

    w8[0] = 0.1012285362903763;
    w8[1] = 0.2223810344533745;
    w8[2] = 0.3137066458778873;
    w8[3] = 0.3626837833783620;
    w8[4] = w8[3];
    w8[5] = w8[2];
    w8[6] = w8[1];
    w8[7] = w8[0];

    g8[0] = -0.9602898564975363;
    g8[1] = -0.7966664774136267;
    g8[2] = -0.5255324099163290;
    g8[3] = -0.1834346424956498;
    g8[4] = -g8[3];
    g8[5] = -g8[2];
    g8[6] = -g8[1];
    g8[7] = -g8[0];

    if (p==1) {x=0.0;w=2.0;}
    if (p==2) {x=g2[m];w=1.0;}
    if (p==3) {x=g3[m];w=w3[m];}
    if (p==4) {x=g4[m];w=w4[m];}
    if (p==5) {x=g5[m];w=w5[m];}  
    if (p==6) {x=g6[m];w=w6[m];} 
    if (p==7) {x=g7[m];w=w7[m];} 
    if (p==8) {x=g8[m];w=w8[m];} 
} // end of line rule function      

// Used by Lobatto 2/3/4D to set Lobatto quadrature points
void line_lobatto_info(
    real_t &x, 
    real_t &w, 
    int &m, 
    int &p){

    m--;
    // creating matrices for weights and points
    real_t L2[2],w3[3],L3[3],w4[4],L4[4],L5[5],w5[5], 
          L6[6],w6[6],L7[7],w7[7],L8[8],w8[8];

    L2[0] = -1.0;
    L2[1] =  1.0;

    w3[0] = 1.0/3.0;
    w3[1] = 4.0/3.0;
    w3[2] = w3[0]; 
        
    L3[0] =  1.0;
    L3[1] =  0.0;
    L3[2] = -1.0; 

    w4[0] = 0.1666666666666666666667;
    w4[1] = 0.8333333333333333333333;
    w4[2] = w4[1]; 
    w4[3] = w4[0]; 

    L4[0] = -1.0;
    L4[1] = -0.447213595499957939282;
    L4[2] = -L4[1]; 
    L4[3] = -L4[0]; 

    L5[0] = -1.0;
    L5[1] = -0.6546536707079771437983;
    L5[2] = 0.0;
    L5[3] = -L5[1]; 
    L5[4] = -L5[0]; 

    w5[0] = 0.1;
    w5[1] = 0.544444444444444444444;
    w5[2] = 0.7111111111111111111111;
    w5[3] = w5[1]; 
    w5[4] = w5[0]; 

    L6[0] = -1.0;
    L6[1] = -0.765055323929464692851;
    L6[2] = -0.2852315164806450963142;
    L6[3] = -L6[2];
    L6[4] = -L6[1];
    L6[5] = -L6[0];


    w6[0] = 0.06666666666666666666667;
    w6[1] = 0.3784749562978469803166;
    w6[2] = 0.5548583770354863530167;
    w6[3] = w6[2];
    w6[4] = w6[1];
    w6[5] = w6[0];

    L7[0] = -1.0;
    L7[1] = -0.830223896278566929872;
    L7[2] = -0.4688487934707142138038;
    L7[3] =  0.0;
    L7[4] = -L7[2];
    L7[5] = -L7[1];
    L7[6] = -L7[0];


    w7[0] = 0.04761904761904761904762;
    w7[1] = 0.276826047361565948011;
    w7[2] = 0.4317453812098626234179;
    w7[3] = 0.487619047619047619048;
    w7[4] = w7[2];
    w7[5] = w7[1];
    w7[6] = w7[0];


    L8[0] = -1.0;
    L8[1] = -0.8717401485096066153375;
    L8[2] = -0.5917001814331423021445;
    L8[3] = -0.2092992179024788687687;
    L8[4] = -L8[3];
    L8[5] = -L8[2];
    L8[6] = -L8[1];
    L8[7] = -L8[0];


    w8[0] = 0.03571428571428571428571;
    w8[1] = 0.210704227143506039383;
    w8[2] = 0.3411226924835043647642;
    w8[3] = 0.4124587946587038815671;
    w8[4] = w8[3];
    w8[5] = w8[2];
    w8[6] = w8[1];
    w8[7] = w8[0];

    if (p==1) {x=0.0;w=2.0;}
    if (p==2) {x=L2[m];w=1.0;}
    if (p==3) {x=L3[m];w=w3[m];}
    if (p==4) {x=L4[m];w=w4[m];}
    if (p==5) {x=L5[m];w=w5[m];}  
    if (p==6) {x=L6[m];w=w6[m];} 
    if (p==7) {x=L7[m];w=w7[m];} 
    if (p==8) {x=L8[m];w=w8[m];} 
} // end of Lobatto line rule function   

// setting gauss quadrature points for 2D elements
void gauss_2d(
    view_c_array <real_t> &these_g_pts,   // gauss points
    view_c_array <real_t> &these_weights, // gauss weights
    view_c_array <real_t> &tot_g_weight,  // 2D product of gauss weights
    int &quad_order){                     // quadrature order (m)

    int tot_pts = (quad_order*quad_order);    // total quad points in 2D

    //find Gauss-Quadrature Points for volume integration
    for (int m = 0; m < tot_pts; m++) {
      
        int i, j;
        real_t x,w;
 
        // sets up the i, j, indices for the line rule implimentation

        j = floor(m/quad_order)+1; 
        i = (m+1) - quad_order*(j-1);

        these_weights(m, 0) = 1.0;
        these_weights(m, 1) = 1.0;
        tot_g_weight(m) = 1.0;

        // xi direction
        line_gauss_info(x,w,i,quad_order); // setting pts/weights in i direction
        these_g_pts(m, 0)  = x; 
        these_weights(m,0)*= w;

        // eta direction
        line_gauss_info(x,w,j,quad_order); // setting pts/weights in j direction
        these_g_pts(m,1)   = x; 
        these_weights(m,1)*= w;

        tot_g_weight(m) = these_weights(m,0)*these_weights(m,1);
    }
} // end Gauss 2D function

// setting gauss quadrature points for 3D elements
void gauss_3d(
    view_c_array <real_t> &these_g_pts,   // gauss points
    view_c_array <real_t> &these_weights, // gauss weights
    view_c_array <real_t> &tot_g_weight,  // 3D product of gauss weights
    int &quad_order){                     // quadrature order (n)

    // total quad points in 3D
    int tot_pts = (quad_order*quad_order*quad_order);  

    //find Gauss-Quadrature Points for volume integration
    for (int m = 0; m < tot_pts; m++) {
      
        int i, j, k, jj;

        real_t x, w;

        // sets up the i, j, k indices for the line rule implimentation
        int p12 = (quad_order * quad_order);

        k=floor(m/p12)+1; 

        jj=(m+1)-p12*(k-1); 
        j=floor((jj-1)/quad_order)+1;

        i=jj-quad_order*(j-1);  

        these_weights(m, 0) = 1.0;
        these_weights(m, 1) = 1.0;
        these_weights(m, 2) = 1.0;
        tot_g_weight(m) = 1.0;

        // xi direction
        line_gauss_info(x,w,i,quad_order); // setting pts/weights in i direction
        these_g_pts(m, 0)   = x; 
        these_weights(m, 0)*= w;

        // eta direction
        line_gauss_info(x,w,j,quad_order); // setting pts/weights in j direction
        these_g_pts(m, 1)   = x; 
        these_weights(m, 1)*= w;

        // mu direction
        line_gauss_info(x,w,k,quad_order); // setting pts/weights in k direction
        these_g_pts(m, 2)   = x; 
        these_weights(m, 2)*= w;

        tot_g_weight(m) = these_weights(m, 0)*these_weights(m, 1)*these_weights(m, 2);
          
    } // end for    
} // end gauss 3D

// setting gauss quadrature points for 4D elements
void gauss_4d(
    view_c_array <real_t> &these_g_pts, // gauss points
    view_c_array <real_t> &these_weights, // gauss weights
    int &quad_order, // quadrature order (n)
    const int &dim){  // dimension

    // total quad points in 3D
    int tot_pts = (quad_order*quad_order*quad_order*quad_order);  

    //find Gauss-Quadrature Points for volume integration
    for (int m = 0; m < tot_pts; m++) {
      
        int i, j, k, l, kk, jj;

        real_t x,w;

        real_t a, b, c;

        a = pow(real_t(quad_order), real_t(dim-1.0));
        b = pow(real_t(quad_order), real_t(dim-2.0));

        // sets up the i, j, k indices for the line rule implimentation
        int p12 = (quad_order * quad_order);

        l = floor(m/(a)) + 1.0;

        kk = (m + 1) - (a)*(l - 1);

        k=floor((kk - 1)/(b)) + 1;

        jj=(m+1) - ((b)*(k - 1)) - (a)*(l - 1);

        j=floor((jj-1)/quad_order)+1;

        i=jj-quad_order*(j-1);  

        these_weights(m, 0) = 1.0;
        these_weights(m, 1) = 1.0;
        these_weights(m, 2) = 1.0;
        these_weights(m, 4) = 1.0;

        // xi direction
        line_gauss_info(x,w,i,quad_order); // setting pts/weights in i direction
        these_g_pts(m, 0)   = x; 
        these_weights(m, 0)*= w;

        // eta direction
        line_gauss_info(x,w,j,quad_order); // setting pts/weights in j direction
        these_g_pts(m, 1)   = x; 
        these_weights(m, 1)*= w;

        // mu direction
        line_gauss_info(x,w,k,quad_order); // setting pts/weights in k direction
        these_g_pts(m, 2)   = x; 
        these_weights(m, 2)*= w;

        // tau direction
        line_gauss_info(x,w,l,quad_order); // setting pts/weights in l direction
        these_g_pts(m, 3)   = x; 
        these_weights(m, 3)*= w;

   } // end for        
} // end Gauss 4D function

// setting Gauss-Lobatto quadrature points for 2D elements
void lobatto_2d(
    view_c_array <real_t> &these_L_pts, // gauss points
    view_c_array <real_t> &these_weights, // gauss weights
    int &quad_order){ // quadrature order (n)

    int tot_pts = (quad_order*quad_order);    // total quad points in 2D

    //find Gauss-Quadrature Points for volume integration
    for (int m = 0; m < tot_pts; m++) {
      
        int i, j;
        real_t x,w;

        // sets up the i, j, indices for the line rule implimentation

        j = floor(m/quad_order)+1; 
        i = (m+1) - quad_order*(j-1);

        these_weights(m, 0) = 1.0;
        these_weights(m, 1) = 1.0;

        // xi direction
        line_lobatto_info(x,w,i,quad_order); // setting pts/weights in i direction
        these_L_pts(m, 0)   = x; 
        these_weights(m, 0)*= w;

        // eta direction
        line_lobatto_info(x,w,j,quad_order); // setting pts/weights in j direction
        these_L_pts(m, 1)   = x; 
        these_weights(m, 1)*= w;
        
    } // end for
} // end function

// setting Gauss-Lobatto quadrature points for 3D elements
void lobatto_3d(
    view_c_array <real_t> &these_L_pts, // gauss points
    view_c_array <real_t> &these_weights, // gauss weights
    int &quad_order){  // quadrature order (n)

    // total quad points in 3D
    int tot_pts = (quad_order*quad_order*quad_order);  

    //find Gauss-Quadrature Points for volume integration
    for (int m = 0; m < tot_pts; m++) {
      
        int i, j, k, jj;

        real_t x,w;

        // sets up the i, j, k indices for the line rule implimentation
        int p12 = (quad_order * quad_order);

        k=floor(m/p12)+1; 

        jj= (m+1)-p12*(k-1); 
        j = floor((jj-1)/quad_order)+1;

        i = jj-quad_order*(j-1);  

        these_weights(m, 0) = 1.0;
        these_weights(m, 1) = 1.0;
        these_weights(m, 2) = 1.0;

        // xi direction
        line_lobatto_info(x,w,i,quad_order); // setting pts/weights in i direction
        these_L_pts(m, 0)   = x; 
        these_weights(m, 0)*= w;

        // eta direction
        line_lobatto_info(x,w,j,quad_order); // setting pts/weights in j direction
        these_L_pts(m, 1)   = x; 
        these_weights(m, 1)*= w;

        // mu direction
        line_lobatto_info(x,w,k,quad_order); // setting pts/weights in k direction
        these_L_pts(m, 2)   = x; 
        these_weights(m, 2)*= w;
          
    } // end for        
} // end function

// setting gauss quadrature points for 4D elements
void lobatto_4d(
    view_c_array <real_t> &these_L_pts, // gauss points
    view_c_array <real_t> &these_weights, // gauss weights
    int &quad_order, // quadrature order (n)
    const int &dim){  // dimension

    // total quad points in 3D
    int tot_pts = (quad_order*quad_order*quad_order*quad_order);  

    //find Gauss-Quadrature Points for volume integration
    for (int m = 0; m < tot_pts; m++) {
      
        int i, j, k, l, kk, jj;

        real_t x,w;

        real_t a, b, c;

        a = pow(real_t(quad_order), real_t(dim-1.0)); //to simplify indexing equations
        b = pow(real_t(quad_order), real_t(dim-2.0));

        // sets up the i, j, k, l indices for the line rule implimentation

        l = floor(m/(a)) + 1.0;

        kk = (m + 1) - (a)*(l - 1);

        k = floor((kk - 1)/(b)) + 1;

        jj = (m+1) - ((b)*(k - 1)) - (a)*(l - 1);

        j = floor((jj-1)/quad_order)+1;

        i = jj-quad_order*(j-1);  

        these_weights(m, 0) = 1.0;
        these_weights(m, 1) = 1.0;
        these_weights(m, 2) = 1.0;
        these_weights(m, 3) = 1.0;

        // xi direction
        line_lobatto_info(x,w,i,quad_order); // setting pts/weights in i direction
        these_L_pts(m, 0)   = x; 
        these_weights(m, 0)*= w;

        // eta direction
        line_lobatto_info(x,w,j,quad_order); // setting pts/weights in j direction
        these_L_pts(m, 1)   = x; 
        these_weights(m, 1)*= w;

        // mu direction
        line_lobatto_info(x,w,k,quad_order); // setting pts/weights in k direction
        these_L_pts(m, 2)   = x; 
        these_weights(m, 2)*= w;

        // tau direction
        line_lobatto_info(x,w,l,quad_order); // setting pts/weights in l direction
        these_L_pts(m, 3)   = x; 
        these_weights(m, 3)*= w;
          
   } // end for        
} // end function

//defining the jacobian for 2D elements
void jacobian_2d(
    view_c_array <real_t> &J_matrix, 
    real_t &det_J,
    const view_c_array <real_t> &vertices, 
    const view_c_array <real_t> &this_partial,
    const int &num_verts){

    int dim = 2;

    // setting jacobian matrix to all zeros
    for(int j = 0; j < dim; j++){  // looping over dimension
        for(int k = 0; k < dim; k++){ //looping over dimension
            J_matrix(j, k) = 0.0;
        }// end for k
    } // end for j

    // solving for the jacobian 
    for(int j = 0; j < dim; j++){  // looping over dimension (partial)
        for(int k = 0; k < dim; k++){ //looping over dimension (vertex position)
            for(int this_x_vert = 0; this_x_vert < num_verts; this_x_vert++){ 

                J_matrix(j, k) += vertices(this_x_vert, k)
                                * this_partial(this_x_vert, j);
                
            } // end for num_verts
        } // end for k
    } // end for j

    //calcualate the determinant of the jacobian for 2D matrices
    det_J = J_matrix(0, 0)*J_matrix(1, 1) 
          - J_matrix(0, 1)*J_matrix(1, 0);
} // end of jacobian_2d function

//defining the jacobian for 3D elements
void jacobian_3d(
    view_c_array <real_t> &J_matrix, 
    real_t &det_J,
    const view_c_array <real_t> &vertices, 
    const view_c_array <real_t> &this_partial,
    const int &num_verts){

    const int dim = 3;

    // setting jacobian matrix to all zeros
    for(int j = 0; j < dim; j++){  // looping over dimension
        for(int k = 0; k < dim; k++){ //looping over dimension
            J_matrix(j, k) = 0.0;
        }// end for k
    } // end for j

    // solving for the jacobian 
    for(int j = 0; j < dim; j++){  // looping over dimension (partial)
        for(int k = 0; k < dim; k++){ //looping over dimension (vertex position)
            //looping over node positions
            for(int vert_lid = 0; vert_lid < num_verts; vert_lid++){ 
               
                J_matrix(j, k) += vertices(vert_lid, k)
                                * this_partial(vert_lid, j);
            
            } // end for num_verts
        } // end for k
    } // end for j



    //calcualate the determinant of the jacobian for 2D/3D matrices
    det_J = J_matrix(0, 0) 
          *(J_matrix(1, 1)*J_matrix(2, 2)) - (J_matrix(2, 1)*J_matrix(1, 2))  //group 1
          - J_matrix(0, 1) 
          *(J_matrix(1, 0)*J_matrix(2, 2)) - (J_matrix(2, 0)*J_matrix(1, 2))  // group 2
          + J_matrix(0, 2) 
          *(J_matrix(1, 0)*J_matrix(2, 1)) - (J_matrix(2, 0)*J_matrix(1, 1)); // group 3 
} // end of jacobian function

//defining the jacobian for 4D elements
void jacobian_4d(
    view_c_array <real_t> &J_matrix, 
    real_t &det_J,
    const view_c_array <real_t> &vertices, 
    const view_c_array <real_t> &this_partial,
    const int &num_verts,
    const int &dim){
   

    // setting jacobian matrix to all zeros
    for(int j = 0; j < dim; j++){  // looping over dimension
        for(int k = 0; k < dim; k++){ //looping over dimension
            J_matrix(j, k) = 0.0;
        }// end for k
    } // end for j
   
    // solving for the jacobian 
    for(int j = 0; j < dim; j++){  // looping over dimension (partial)
        for(int k = 0; k < dim; k++){ //looping over dimension (vertex position)
        //looping over node positions
            for(int this_x_vert = 0; this_x_vert < num_verts; this_x_vert++){ 
                J_matrix(j, k) += vertices(this_x_vert, k)
                                * this_partial(this_x_vert, j);
            } // end for num_verts
        } // end for k
    } // end for j
    
    det_J = J_matrix(0, 3) * J_matrix(1, 2) * J_matrix(2, 1) * J_matrix(3, 0) 
          - J_matrix(0, 2) * J_matrix(1, 3) * J_matrix(2, 1) * J_matrix(3, 0) 
          - J_matrix(0, 3) * J_matrix(1, 1) * J_matrix(2, 2) * J_matrix(3, 0) 
          + J_matrix(0, 1) * J_matrix(1, 3) * J_matrix(2, 2) * J_matrix(3, 0) 
          + J_matrix(0, 2) * J_matrix(1, 1) * J_matrix(2, 3) * J_matrix(3, 0) 
          - J_matrix(0, 1) * J_matrix(1, 2) * J_matrix(2, 3) * J_matrix(3, 0) 
          - J_matrix(0, 3) * J_matrix(1, 2) * J_matrix(2, 0) * J_matrix(3, 1) 
          + J_matrix(0, 2) * J_matrix(1, 3) * J_matrix(2, 0) * J_matrix(3, 1) 
          + J_matrix(0, 3) * J_matrix(1, 0) * J_matrix(2, 2) * J_matrix(3, 1) 
          - J_matrix(0, 0) * J_matrix(1, 3) * J_matrix(2, 2) * J_matrix(3, 1) 
          - J_matrix(0, 2) * J_matrix(1, 0) * J_matrix(2, 3) * J_matrix(3, 1) 
          + J_matrix(0, 0) * J_matrix(1, 2) * J_matrix(2, 3) * J_matrix(3, 1) 
          + J_matrix(0, 3) * J_matrix(1, 1) * J_matrix(2, 0) * J_matrix(3, 2) 
          - J_matrix(0, 1) * J_matrix(1, 3) * J_matrix(2, 0) * J_matrix(3, 2) 
          - J_matrix(0, 3) * J_matrix(1, 0) * J_matrix(2, 1) * J_matrix(3, 2) 
          + J_matrix(0, 0) * J_matrix(1, 3) * J_matrix(2, 1) * J_matrix(3, 2) 
          + J_matrix(0, 1) * J_matrix(1, 0) * J_matrix(2, 3) * J_matrix(3, 2) 
          - J_matrix(0, 0) * J_matrix(1, 1) * J_matrix(2, 3) * J_matrix(3, 2) 
          - J_matrix(0, 2) * J_matrix(1, 1) * J_matrix(2, 0) * J_matrix(3, 3) 
          + J_matrix(0, 1) * J_matrix(1, 2) * J_matrix(2, 0) * J_matrix(3, 3) 
          + J_matrix(0, 2) * J_matrix(1, 0) * J_matrix(2, 1) * J_matrix(3, 3) 
          - J_matrix(0, 0) * J_matrix(1, 2) * J_matrix(2, 1) * J_matrix(3, 3) 
          - J_matrix(0, 1) * J_matrix(1, 0) * J_matrix(2, 2) * J_matrix(3, 3) 
          + J_matrix(0, 0) * J_matrix(1, 1) * J_matrix(2, 2) * J_matrix(3, 3);
} // end of jacobian function

//defining the inverse jacobian for 2D element    
void jacobian_inverse_2d(
    view_c_array <real_t> &J_inverse, 
    const view_c_array <real_t> &jacobian){

    real_t det = 0.0;
    det = jacobian(0, 0)*jacobian(1, 1) 
        - jacobian(0, 1)*jacobian(1, 0);


    J_inverse(0, 0) =  jacobian(1, 1)/det;
    J_inverse(0, 1) = -jacobian(0, 1)/det;
    J_inverse(1, 0) = -jacobian(1, 0)/det;
    J_inverse(1, 1) =  jacobian(0, 0)/det;
} // end of 2D jacobin inverse

// defining  the inverse of the Jacobian for 3D elements
void jacobian_inverse_3d(
    view_c_array <real_t> &J_inverse,
    const view_c_array <real_t> &jacobian){

    real_t A_11 = jacobian(1, 1)*jacobian(2, 2) 
                - jacobian(1, 2)*jacobian(2, 1);

    real_t A_22 = jacobian(2, 2)*jacobian(0, 0) 
                - jacobian(2, 0)*jacobian(0, 2);

    real_t A_33 = jacobian(0, 0)*jacobian(1, 1) 
                - jacobian(0, 1)*jacobian(1, 0);

    real_t A_12 = jacobian(1, 2)*jacobian(2, 0) 
                - jacobian(1, 0)*jacobian(2, 2);

    real_t A_23 = jacobian(2, 0)*jacobian(0, 1) 
                - jacobian(2, 1)*jacobian(0, 0);

    real_t A_31 = jacobian(0, 1)*jacobian(1, 2) 
                - jacobian(0, 2)*jacobian(1, 1);

    real_t A_21 = jacobian(2, 1)*jacobian(0, 2) 
                - jacobian(0, 1)*jacobian(2, 2);

    real_t A_32 = jacobian(0, 2)*jacobian(1, 0) 
                - jacobian(1, 2)*jacobian(0, 0);

    real_t A_13 = jacobian(1, 0)*jacobian(1, 1) 
                - jacobian(2, 0)*jacobian(1, 1);

    real_t  det = jacobian(0, 0)*A_11 + jacobian(0, 1)*A_21 
                + jacobian(0, 2)*A_31;


    J_inverse(0, 0) = A_11/det;
    J_inverse(0, 1) = A_12/det;
    J_inverse(0, 2) = A_13/det;
    J_inverse(1, 0) = A_21/det;
    J_inverse(1, 1) = A_22/det;
    J_inverse(1, 2) = A_23/det;
    J_inverse(2, 0) = A_31/det;
    J_inverse(2, 1) = A_32/det;
    J_inverse(2, 2) = A_33/det;
} // end of inverse jacobian 

// defining  the inverse of the Jacobian for 4D elements
void jacobian_inverse_4d(
    view_c_array <real_t> &J_inverse,
    const view_c_array <real_t> &jacobian,
    const real_t &det_J){


    real_t inv[16];
    real_t m[16];
    real_t inv_mat[16];

    // convering jacobian into array length 16
    for (int array = 0; array < 4; array++){ 
        m[array] = jacobian(0, array);
    }
    for (int array = 4; array < 8; array++){ 
        m[array] = jacobian(1, array-4);
    }
    for (int array = 8; array < 12; array++){ 
        m[array] = jacobian(2, array-8);
    }
    for (int array = 12; array < 16; array++){ 
        m[array] = jacobian(3, array-12);     
    }


    inv[0] = m[5]  * m[10] * m[15] 
           - m[5]  * m[11] * m[14] 
           - m[9]  * m[6]  * m[15] 
           + m[9]  * m[7]  * m[14] 
           + m[13] * m[6]  * m[11] 
           - m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] 
           + m[4]  * m[11] * m[14] 
           + m[8]  * m[6]  * m[15] 
           - m[8]  * m[7]  * m[14] 
           - m[12] * m[6]  * m[11] 
           + m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] 
           - m[4]  * m[11] * m[13] 
           - m[8]  * m[5] * m[15] 
           + m[8]  * m[7] * m[13] 
           + m[12] * m[5] * m[11] 
           - m[12] * m[7] * m[9];

    inv[12]=-m[4]  * m[9] * m[14] 
           + m[4]  * m[10] * m[13] 
           + m[8]  * m[5] * m[14] 
           - m[8]  * m[6] * m[13] 
           - m[12] * m[5] * m[10] 
           + m[12] * m[6] * m[9];

    inv[1] =-m[1]  * m[10] * m[15] 
           + m[1]  * m[11] * m[14] 
           + m[9]  * m[2] * m[15] 
           - m[9]  * m[3] * m[14] 
           - m[13] * m[2] * m[11] 
           + m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] 
           - m[0]  * m[11] * m[14] 
           - m[8]  * m[2] * m[15] 
           + m[8]  * m[3] * m[14] 
           + m[12] * m[2] * m[11] 
           - m[12] * m[3] * m[10];

    inv[9] =-m[0]  * m[9] * m[15] 
           + m[0]  * m[11] * m[13] 
           + m[8]  * m[1] * m[15] 
           - m[8]  * m[3] * m[13] 
           - m[12] * m[1] * m[11] 
           + m[12] * m[3] * m[9];

    inv[13]= m[0]  * m[9] * m[14] 
           - m[0]  * m[10] * m[13] 
           - m[8]  * m[1] * m[14] 
           + m[8]  * m[2] * m[13] 
           + m[12] * m[1] * m[10] 
           - m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] 
           - m[1]  * m[7] * m[14] 
           - m[5]  * m[2] * m[15] 
           + m[5]  * m[3] * m[14] 
           + m[13] * m[2] * m[7] 
           - m[13] * m[3] * m[6];

    inv[6] =-m[0]  * m[6] * m[15] 
           + m[0]  * m[7] * m[14] 
           + m[4]  * m[2] * m[15] 
           - m[4]  * m[3] * m[14] 
           - m[12] * m[2] * m[7] 
           + m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] 
           - m[0]  * m[7] * m[13] 
           - m[4]  * m[1] * m[15] 
           + m[4]  * m[3] * m[13] 
           + m[12] * m[1] * m[7] 
           - m[12] * m[3] * m[5];

    inv[14]=-m[0]  * m[5] * m[14] 
           + m[0]  * m[6] * m[13] 
           + m[4]  * m[1] * m[14] 
           - m[4]  * m[2] * m[13] 
           - m[12] * m[1] * m[6] 
           + m[12] * m[2] * m[5];

    inv[3] =-m[1] * m[6] * m[11] 
           + m[1] * m[7] * m[10] 
           + m[5] * m[2] * m[11] 
           - m[5] * m[3] * m[10] 
           - m[9] * m[2] * m[7] 
           + m[9] * m[3] * m[6];

    inv[7] =  m[0] * m[6] * m[11] 
           - m[0] * m[7] * m[10] 
           - m[4] * m[2] * m[11] 
           + m[4] * m[3] * m[10] 
           + m[8] * m[2] * m[7] 
           - m[8] * m[3] * m[6];

    inv[11]= -m[0] * m[5] * m[11] 
           +  m[0] * m[7] * m[9] 
           +  m[4] * m[1] * m[11] 
           -  m[4] * m[3] * m[9] 
           -  m[8] * m[1] * m[7] 
           +  m[8] * m[3] * m[5];

    inv[15]= m[0] * m[5] * m[10] 
           - m[0] * m[6] * m[9] 
           - m[4] * m[1] * m[10] 
           + m[4] * m[2] * m[9] 
           + m[8] * m[1] * m[6] 
           - m[8] * m[2] * m[5];

    real_t det_inv = (1.0 / det_J);

    for (int i = 0; i < 16; i++){
        inv_mat[i] = inv[i] * det_inv;
    }
   
   // convering inverse jacobian back to 2D array
    for (int array = 0; array < 4; array++){ 
        J_inverse(0, array) = inv_mat[array];
    }

    for (int array = 4; array < 8; array++){
        J_inverse(1, array-4) = inv_mat[array];      
    }

    for (int array = 8; array < 12; array++){
        J_inverse(2, array-8) = inv_mat[array];    
    }

    for (int array = 12; array < 16; array++){ 
        J_inverse(3, array-12) = inv_mat[array];          
    }
} // end of quad point inverse jacobian functions

// creates nodal positions with Chebyshev spacing
void chebyshev_nodes_1D(
    view_c_array <real_t> &cheb_nodes_1D,   // Chebyshev nodes
    const int &order){                      // Interpolation order

    real_t pi = 3.14159265358979323846;

    // loop to build the base set of Chebyshev nodes
    for (int i = 1; i < order + 1; i++){
        cheb_nodes_1D(i - 1) = -cos(pi*(2.0*i - 1.0)/(2.0*(order + 1.0)));
    } // end for i

    //replacing the first and last index with the element boundary 
    cheb_nodes_1D(0) = -1.0;
    cheb_nodes_1D(order) = 1.0;
}



/*
 .-------------------------------. 
| .----------------------------. |
| |    _____       ________    | |
| |   / ___ `.    |_   ___ `.  | |
| |  |_/___) |      | |   `. \ | |
| |   .'____.'      | |    | | | |
| |  / /____       _| |___.' / | |
| |  |_______|    |________.'  | |
| |                            | |
| '----------------------------' |
 '-------------------------------' 
*/

/*
===========================
2D Quad 4 Elements
===========================


The finite element local point numbering for a 4 node Hexahedral is
as follows

        Eta
         ^
         |
  3------+-----2
  |      |     |
  |      |     |
  |      |     |
  |      ------+------> Xi   
  |            |
  |            |
  0------------1

*/

real_t Quad4::ref_vert[Quad4::num_verts*Quad4::num_dim] = 
    
    {// listed as {Xi, Eta}
    -1.0, -1.0,// 0
     1.0, -1.0,// 1
     1.0,  1.0,// 2
    -1.0,  1.0,// 3
    };

const int Quad4::vert_to_node[Quad4::num_verts] = 
    {
    0,
    2,
    6,
    8
    };

// calculate a physical position in an element for a given xi,eta
void Quad4::physical_position(
    view_c_array <real_t> &x_point, 
    const view_c_array <real_t> &xi_point, 
    const view_c_array <real_t> &vertices){

    real_t basis_a[num_verts];
    auto basis = view_c_array <real_t> (basis_a, num_verts);

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the shape functions from each vertex for 0 through num_verts(xi,eta)
    for( int vert_lid = 0; vert_lid < num_verts; vert_lid++ ){
        basis(vert_lid) = 1.0/4.0
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1));
    }// end for

    // calculate the position in physical space
    for (int dim = 0; dim < num_dim; dim++) x_point(dim) = 0.0;
    
    for (int vert_lid = 0; vert_lid < num_verts; vert_lid++ ){
        for (int dim = 0; dim < num_dim; dim++){
            x_point(dim) += vertices(vert_lid, dim)*basis(vert_lid);
        }// end for dim
    } // end for vert_lid

} // end of physical position functionfunction


// calculate the value for the basis at each node for a given xi,eta
void Quad4::basis(
    view_c_array <real_t> &basis,
    const view_c_array <real_t> &xi_point){

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);
    
    // calculate the shape functions from each vertex for 0 through num_verts(xi,eta)
    for( int vert_lid = 0; vert_lid < num_verts; vert_lid++ ){

        basis(vert_lid) = (1.0/4.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1));

    }// end for

}// end of quad4 basis functions 


// Partial derivative of shape functions with respect to Xi
void  Quad4::partial_xi_basis(
    view_c_array <real_t>  &partial_xi, 
    const view_c_array <real_t> &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);
    
    for (int vert_lid = 0; vert_lid < num_verts; vert_lid++){
        partial_xi(vert_lid) = (1.0/4.0)
            * (ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1));
    }

}// end of partial xi funciton


// Partial derivative of shape functions with respect to Eta
void  Quad4::partial_eta_basis(
    view_c_array <real_t> &partial_eta, 
    const view_c_array <real_t> &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    for (int vert_lid = 0; vert_lid < num_verts; vert_lid++){
        partial_eta(vert_lid) = (1.0/4.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (ref_verts(vert_lid, 1));
    }

}// end of partial eta function

inline int Quad4::vert_node_map(const int vert_lid){

    return vert_to_node[vert_lid];

}



/*
===========================
2D Quad 8 Elements
===========================


 The finite element local point numbering for a 8 node Hexahedral is
 as follows

         Eta
          ^
          |
  3-------6------2
  |       |      |
  |       |      |
  |       |      |
  |       |      |
  7       +------5-----> Xi   
  |              |
  |              |
  |              |
  0------4-------1

*/

real_t Quad8::ref_vert[Quad8::num_verts*Quad8::num_dim] = // listed as {Xi, Eta}
    {// listed as {Xi, Eta}
    -1.0, -1.0, // 0  
     1.0, -1.0, // 1
     1.0,  1.0, // 2
    -1.0,  1.0, // 3
    // midline nodes
     0.0, -1.0, // 4
     1.0,  0.0, // 5
     0.0,  1.0, // 6
    -1.0,  0.0, // 7
    };

const int Quad8::vert_to_node[Quad8::num_verts] = 
    {
    0,
    4,
    24,
    20,
    2,
    14,
    23,
    10
    };

// calculate a physical position in an element for a given xi,eta,
void Quad8::physical_position(
    view_c_array <real_t> &x_point, 
    const view_c_array <real_t> &xi_point, 
    const view_c_array <real_t> &vertices){

    real_t basis_a[num_verts];
    auto basis = view_c_array <real_t> (basis_a, num_verts);

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the shape functions for node 0,1,2,3(xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
        basis(vert_lid) = (1.0/4.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (xi_point(0)*ref_verts(vert_lid, 0) 
            +  xi_point(1)*ref_verts(vert_lid, 1) - 1.0);
    } // end for vert_lid

    // calculate the shape functions for node 4,6(xi,eta)
    for( int vert_lid = 4; vert_lid <= 6; vert_lid += 2 ){
        basis(vert_lid) = (1.0/2.0)
            * (1.0 - xi_point(0)*xi_point(0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1));
    } // end for vert_lid

    // calculate the shape functions for node 5,7 (xi,eta)
    for( int vert_lid = 5; vert_lid <= 7; vert_lid += 2 ){
        basis(vert_lid) = (1.0/2.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 - xi_point(1)*xi_point(1));
    } // end for vert_lid

    // calculate the position in physical space
    for (int dim = 0; dim < num_dim; dim++) x_point(dim) = 0.0;

    for (int vert_lid = 0; vert_lid < num_verts; vert_lid++ ){
        for (int dim = 0; dim < num_dim; dim++){
            x_point(dim) += vertices(vert_lid, dim)*basis(vert_lid);
        } // end for dim
    } // end for vert_lid

} // end of function


// calculate the value for the basis at each node for a given xi,eta
void Quad8::basis(
    view_c_array <real_t> &basis,
    const view_c_array <real_t> &xi_point){

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);
    
    // calculate the shape functions for node 0,1,2,3(xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
        basis(vert_lid) = (1.0/4.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (xi_point(0)*ref_verts(vert_lid, 0) 
            +  xi_point(1)*ref_verts(vert_lid, 1) - 1.0);
    } // end for vert_lid


    // calculate the shape functions for node 4,6(xi,eta)
    for( int vert_lid = 4; vert_lid <= 6; vert_lid += 2 ){
        basis(vert_lid) = (1.0/2.0)
            * (1.0 - xi_point(0)*xi_point(0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1));
    } // end for vert_lid

    // calculate the shape functions for node 5,7 (xi,eta)
    for( int vert_lid = 5; vert_lid <= 7; vert_lid += 2 ){
        basis(vert_lid) = (1.0/2.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 - xi_point(1)*xi_point(1));

    } // end for vert_lid

}// end of quad8 basis functions


// Partial derivative of shape functions with respect to Xi
void Quad8::partial_xi_basis(
    view_c_array <real_t>  &partial_xi, 
    const view_c_array <real_t> &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the Xi partials for node 0,1,2,3 (xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
        partial_xi(vert_lid) = 1.0/4.0
            * (ref_verts(vert_lid, 0))
            * (1.0 + ref_verts(vert_lid, 1)*xi_point(1))
            *((2.0 * ref_verts(vert_lid, 0)*xi_point(0)) 
            + (ref_verts(vert_lid, 1)*xi_point(1)));
    } // end for vert_lid


    // calculate the Xi partials for node 4,6 (xi,eta)
    for( int vert_lid = 4; vert_lid <= 6; vert_lid += 2 ){
        partial_xi(vert_lid) = -1.0
            * (xi_point(0))
            * (1 + ref_verts(vert_lid, 1)*xi_point(1));
    } // end for vert_lid

    // calculate the Xi partials for node 5,7 (xi,eta)
    for( int vert_lid = 5; vert_lid <= 7; vert_lid += 2 ){
        partial_xi(vert_lid) = 1.0/2.0
            * (ref_verts(vert_lid, 0))
            * (1.0 - xi_point(1)*xi_point(1));

    } // end for vert_lid

} // end partial Xi function


// Partial derivative of shape functions with respect to Eta
void Quad8::partial_eta_basis(
    view_c_array <real_t>  &partial_eta, 
    const view_c_array <real_t> &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the Eta partials for node 0,1,2,3 (xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
        partial_eta(vert_lid) = (1.0/4.0)
            * (1.0 + ref_verts(vert_lid, 0)*xi_point(0))
            * (ref_verts(vert_lid, 1))
            *((ref_verts(vert_lid, 0)*xi_point(0))
            + (2.0 * ref_verts(vert_lid, 1)*xi_point(1))); 
    } // end for vert_lid

    // calculate the Eta partials for node 4,6 (xi,eta)
    for( int vert_lid = 4; vert_lid <= 6; vert_lid += 2 ){
        partial_eta(vert_lid) = (1.0/2.0)
            * (1.0 - xi_point(0)*xi_point(0))
            * (ref_verts(vert_lid, 1));
   } // end for vert_lid

    // calculate the Eta partials for node 5,7 (xi,eta)
    for( int vert_lid = 5; vert_lid <= 7; vert_lid += 2 ){
        partial_eta(vert_lid) = (-1.0)
            * (1.0 + ref_verts(vert_lid, 0)*xi_point(0))
            * (xi_point(1));
    } // end for vert_lid

} // end partial Eta function


inline int Quad8::vert_node_map(const int vert_lid){

    return vert_to_node[vert_lid];

}




/*
===========================
2D Quad 12 Elements
===========================


 The finite element local point numbering for a 8 node Hexahedral is
 as follows (NEED TO DEFINE)

         Eta
          ^
          |
  3---7------6---2
  |       |      |
  |       |      |
 11       |      10
  |       |      |
  |       +------|-----> Xi   
  |              |
  8              9
  |              |
  0----4-----5---1

*/
      
real_t Quad12::ref_vert[Quad12::num_verts*Quad12::num_dim] = 
    {// listed as {Xi, Eta}
    //corner nodes
    -1.0, -1.0 ,// 0
     1.0, -1.0 ,// 1
     1.0,  1.0 ,// 2
    -1.0,  1.0 ,// 3
    // Eta +- 1./3.
    -1./3., -1.0 ,// 4
     1./3., -1.0 ,// 5
     1./3.,  1.0 ,// 6
    -1./3.,  1.0 ,// 7
    // Xi +- 1./3.
    -1.0, -1./3. ,// 8
     1.0, -1./3. ,// 9
     1.0,  1./3. ,// 10
    -1.0,  1./3. ,// 11
    };

const int Quad12::vert_to_node[Quad12::num_verts] = 
    {
    0,
    6,
    48,
    42,
    2,
    4,
    46,
    44,
    14,
    20,
    34,
    28
    };

// calculate a physical position in an element for a given xi,eta,
void Quad12::physical_position(
    view_c_array <real_t>  &x_point, 
    const view_c_array <real_t>  &xi_point, 
    const view_c_array <real_t>  &vertices){

    real_t basis_a[num_verts];
    auto basis = view_c_array <real_t> (basis_a, num_verts);


    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the shape functions for node 0,1,2,3(xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
        basis(vert_lid) = 1.0/32.0
            * (1.0 + ref_verts(vert_lid, 0) * xi_point(0))
            * (1.0 + ref_verts(vert_lid, 1) * xi_point(1))
            * (9.0 * (xi_point(0) * xi_point(0) + xi_point(1) * xi_point(1))
            - (10.0));

    } // end for vert_lid

    // calculate the shape functions for node 4-7(xi,eta)
    for( int vert_lid = 4; vert_lid <= 7; vert_lid++ ){
        basis(vert_lid) = 9.0/32.0
            * (1.0 - xi_point(0) * xi_point(0))
            * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
            * (1.0 + 9.0 * xi_point(0) * ref_verts(vert_lid, 0));
    } // end for vert_lid

    // calculate the shape functions for node 8-11 (xi,eta)
    for( int vert_lid = 8; vert_lid <= 11; vert_lid++ ){
        basis(vert_lid) = 9.0/32.0
            * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
            * (1.0 - xi_point(1) * xi_point(1))
            * (1.0 + 9.0 * xi_point(1) * ref_verts(vert_lid, 1));
    } // end for vert_lid

    // calculate the position in physical space
    for (int dim = 0; dim < num_dim; dim++){
        x_point(dim) = 0.0;
    }

    for (int vert_lid = 0; vert_lid < num_verts; vert_lid++ ){
        for (int dim = 0; dim < num_dim; dim++){
            x_point(dim) += vertices(vert_lid, dim)*basis(vert_lid);
        } // end for dim
    } // end for vert_lid

} // end of function


// calculate the value for the basis at each node for a given xi,eta
void Quad12::basis(
    view_c_array <real_t>  &basis,
    const view_c_array <real_t>  &xi_point){

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the shape functions for node 0,1,2,3(xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
        basis(vert_lid) = 1.0/32.0
            * (1.0 + ref_verts(vert_lid, 0) * xi_point(0))
            * (1.0 + ref_verts(vert_lid, 1) * xi_point(1))
            * (9.0 * (xi_point(0) * xi_point(0) + xi_point(1) * xi_point(1))
            - (10.0));

    } // end for vert_lid

    // calculate the shape functions for node 4-7(xi,eta)
    for( int vert_lid = 4; vert_lid <= 7; vert_lid++ ){
        basis(vert_lid) = 9.0/32.0
            * (1.0 - xi_point(0) * xi_point(0))
            * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
            * (1.0 + 9.0 * xi_point(0) * ref_verts(vert_lid, 0));
    } // end for vert_lid

    // calculate the shape functions for node 8-11 (xi,eta)
    for( int vert_lid = 8; vert_lid <= 11; vert_lid++ ){
        basis(vert_lid) = 9.0/32.0
                         * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
                         * (1.0 - xi_point(1) * xi_point(1))
                         * (1.0 + 9.0 * xi_point(1) * ref_verts(vert_lid, 1));
    } // end for vert_lid

}// end of quad12 basis functions


// Partial derivative of shape functions with respect to Xi
void Quad12::partial_xi_basis(
    view_c_array <real_t>  &partial_xi, 
    const view_c_array <real_t>  &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the Xi partials for node 0,1,2,3 (xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
        partial_xi(vert_lid) = 1.0/32.0
            * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
            *((9.0 * ref_verts(vert_lid, 0) 
            * (xi_point(0) * xi_point(0) + xi_point(1) * xi_point(1)))
            + (18.0 * xi_point(0) * (1.0 + xi_point(0) * ref_verts(vert_lid, 0)))
            - (10.0 * ref_verts(vert_lid, 0)));
    } // end for vert_lid

    // calculate the Xi partials for node 4,5,6,7 (xi,eta)
    for( int vert_lid = 4; vert_lid < 8; vert_lid++ ){
        partial_xi(vert_lid) = (9.0/32.0) 
                              * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
                              *((9.0 * ref_verts(vert_lid, 0) 
                              * (1.0 - 3.0 * xi_point(0)*xi_point(0)))
                              - (2.0 * xi_point(0)));
    } // end for vert_lid

    // calculate the Xi partials for node 8,9,10,11 (xi,eta)
    for( int vert_lid = 8; vert_lid <= 11; vert_lid++){
        partial_xi(vert_lid) = 9.0/32.0
            * (ref_verts(vert_lid, 0))
            * (1.0 - xi_point(1) * xi_point(1))
            * (1.0 + 9.0 * xi_point(1) * ref_verts(vert_lid, 1));
    } // end for vert_lid

} // end partial Xi function


// Partial derivative of shape functions with respect to Eta
void Quad12::partial_eta_basis(
    view_c_array <real_t> &partial_eta, 
    const view_c_array <real_t>  &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);
    // calculate the Eta partials for node 0,1,2,3 (xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
        partial_eta(vert_lid) = 1.0/32.0
            * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
            *((9.0 * ref_verts(vert_lid, 1) 
            * (xi_point(0) * xi_point(0) + xi_point(1) * xi_point(1)))
            + (18.0 * xi_point(1) * (1.0 + xi_point(1) * ref_verts(vert_lid, 1)))
            - (10.0 * ref_verts(vert_lid, 1)));
    } // end for vert_lid

    // calculate the Eta partials for node 4,5,6,7 (xi,eta)
    for( int vert_lid = 4; vert_lid <= 7; vert_lid++ ){
        partial_eta(vert_lid) = 9.0/32.0
            * (1.0 - xi_point(0) * xi_point(0))
            * (1.0 + 9.0 * xi_point(0) * ref_verts(vert_lid, 0))
            * (ref_verts(vert_lid, 1));
    } // end for vert_lid

    // calculate the Eta partials for node 8,9,10,11 (xi,eta)
    for( int vert_lid = 8; vert_lid <= 11; vert_lid++){
        partial_eta(vert_lid) = 9.0/32.0
            * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
            *((9.0 * ref_verts(vert_lid, 1) * (1.0 - 3.0 * xi_point(1)*xi_point(1)))
            - (2.0 * xi_point(1)));

    } // end for vert_lid

} // end partial Eta function

inline int Quad12::vert_node_map(const int vert_lid){

    return vert_to_node[vert_lid];
}

/*
===========================
2D Quad 12 Elements
===========================


 The finite element local point numbering for a 8 node Hexahedral is
 as follows (NEED TO DEFINE)

         Eta
          ^
          |
  3---7------6---2
  |       |      |
  |       |      |
 11       |      10
  |       |      |
  |       +------|-----> Xi   
  |              |
  8              9
  |              |
  0----4-----5---1

*/



/*
 ==========================
  Arbitrary Order Elements
 ==========================


  / __ \                | | \ | |
 | |  | |_   _  __ _  __| |  \| |
 | |  | | | | |/ _` |/ _` | . ` |
 | |__| | |_| | (_| | (_| | |\  |
  \___\_\\__,_|\__,_|\__,_|_| \_| 

representative linear element for visualization
   
         Eta (j)
          ^
          |
  3--------------2
  |       |      |
  |       |      |
  |       |      |
  |       |      |
  |       +------|-----> Xi (i) 
  |              |
  |              |
  |              |
  0--------------1
*/


// Lagrange Interp in 1D, returns interpolants and derivative
// works with any nodal spacing
void QuadN::lagrange_1D(
    view_c_array <real_t> &interp,          // interpolant
    view_c_array <real_t> &Dinterp,         // derivative of function
    const real_t &x_point,                  // point of interest in element
    const view_c_array <real_t> &xi_point,  // nodal positions in 1D, normally chebyshev
    const int &orderN){                     // order of element

    real_t num_a[orderN+1];
    auto num = view_c_array <real_t> (num_a, orderN+1); // numerator of interpolant

    real_t denom_a[orderN+1];
    auto denom = view_c_array <real_t> (denom_a, orderN+1); // denomenator of interpolant
  
    real_t q = 0.0;
   
    for(int i = 0; i < orderN + 1; i++){ // looping over the nodes
        real_t n = 1.0;         // placeholder numerator
        real_t d = 1.0;         // placeholder denominator
        real_t c = 1.0;         // placeholder value of n/d
        real_t p = 0.0;         // placeholder for derivative values
        real_t s = 1.0;
        
        for(int j = 0; j < orderN + 1; j++){  // looping over the nodes !=i
            if (j != i ){
                n = n*(x_point - xi_point(j));
                d = d*(xi_point(i) - xi_point(j));
                real_t s = 1.0;
                
                for(int N = 0; N < orderN + 1; N++){  // looping over the nodes !=i
                    if (N != j && N != i ){
                    s = s * (x_point - xi_point(N));
                    }// end if
                }//end for
                
                p += s; 
            }//end if
         
            c = n/d; // storing a single value for interpolation for node i
            q = (p/d); // storing the derivative of the interpolating function
        } // end looping over nodes != i

        // writing value to vectors for later use
        interp(i)   = c;     // Interpolant value at given point
        Dinterp(i)  = q;     // derivative of each function
    } // end loop over all nodes

} // end of Legrange_1D function


// Corners of Lagrange element for mapping
void QuadN::corners (
    view_c_array <real_t> &lag_nodes,   // Nodes of Lagrange elements 
    view_c_array <real_t> &lag_corner,  // corner nodes of HexN element
    const int &orderN){                 // Element order

    /*
    This image represents the corner mapping notation of an arbitrary ordered
    Lagrange element. The corner function takes in the element order and nodal positions and
    returns a vector containing the indices of the corner in alphabetical order.

          Eta
           ^
           |
    C------+-----D
    |      |     |
    |      |     |
    |      |     |
    |      ------+------> Xi   
    |            |
    |            |
    A------------B

    */

    int num_corners = 4;
    int N = orderN + 1;      //number of nodes in each direction
    int corner_ids[num_corners];

    corner_ids[0] = 0;                      
    corner_ids[1] = N - 1.0;                 
    corner_ids[2] = (N*N) - N;         
    corner_ids[3] = (N*N)-1.0;               

    for(int corner = 0; corner < num_corners; corner++){
        for(int dim = 0; dim < num_dim; dim++){
            lag_corner(corner, dim) = lag_nodes(corner_ids[corner], dim);
        }
    }

    // lag_corner[0] = lag_nodes[A];
    // std::cout<<"Corner A = "<<lag_corner[0][0] << "  "<<lag_corner[0][1] <<std::endl;

    // lag_corner[1] = lag_nodes[B];
    // std::cout<<"Corner B = "<<lag_corner[1][0] << "  "<<lag_corner[1][1] <<std::endl;

    // lag_corner[2] = lag_nodes[C];
    // std::cout<<"Corner C = "<<lag_corner[2][0] << "  "<<lag_corner[2][1] <<std::endl;

    // lag_corner[3] = lag_nodes[D];
    // std::cout<<"Corner D = "<<lag_corner[3][0] << "  "<<lag_corner[3][1] <<std::endl;

}// end of corner mapping function


// Functions for mapping reference position to physical position for any 
// point in an arbitrary order 3D lagrange element
void QuadN::physical_position (
    view_c_array <real_t> &x_point,             // location in real space
    const view_c_array <real_t> &lag_nodes,     // Nodes of Lagrange elements 
    const view_c_array <real_t> &lag_basis_2d,  // 2D basis values 
    const int &orderN){                         // order of the element

    int nodes = orderN + 1;
    int Nnodes_2d = nodes * nodes;

    for (int this_vert = 0; this_vert < Nnodes_2d; this_vert++ ){
        for (int dim = 0; dim < num_dim; dim++){
            x_point(dim) += lag_nodes(this_vert, dim)*lag_basis_2d(this_vert);
        } // end for dim
    } // end for this_vert

}// end physical position function


void QuadN::basis_partials (
    view_c_array <real_t> &lag_nodes,       // Nodes of Lagrange elements (to be filled in)
    view_c_array <real_t> &nodes_1d,        // Nodal spacing in 1D, any spacing is accepted
    view_c_array <real_t> &val_1d,          // Interpolant Value in 1D
    view_c_array <real_t> &DVal_1d,         // Derivateive of basis in 1D
    view_c_array <real_t> &val_2d,          // for holding the interpolant in each direction
    view_c_array <real_t> &DVal_2d,         // for holding the derivatives in each direction
    view_c_array <real_t> &lag_basis_2d,    // 3D basis values 
    view_c_array <real_t> &lag_partial,     // Partial of basis 
    const view_c_array <real_t> &xi_point,  // point of interest
    const int &orderN){                     // Element order

    /*

    representative linear element for visualization

          Eta
           ^
           |
    3------+-----2
    |      |     |
    |      |     |
    |      |     |
    |      ------+------> Xi   
    |            |
    |            |
    0------------1


    */

    int N = orderN + 1;      //number of nodes in each direction
    int tot_pts = (N*N);     // total nodes in 2D

    real_t sumi = 0.0;
    real_t sumj = 0.0;

    //Setting nodes for Lagrange Elements
    for (int m = 0; m < tot_pts; m++) {

        int i, j;
        // sets up the i and j indices for the nodes of an 
        j = floor(m/N)+1; 
        i = (m+1) - N*(j-1);

        // xi direction
        lag_nodes(m, 0) = nodes_1d(i-1); 

        // eta direction
        lag_nodes(m, 1) = nodes_1d(j-1); 


        // calling function to assign nodal values for basis and derivative
        
        //evaluating Lagrange interpolants for each function at xi_point
        lagrange_1D(val_1d, DVal_1d, xi_point(0), nodes_1d, orderN);
        val_2d(m, 0)  = val_1d(i-1); 
        DVal_2d(m, 0) = DVal_1d(i-1);

        // resetting to zero
        for(int i = 0.0; i < N; i++){
            val_1d(i)  = 0.0;
            DVal_1d(i) = 0.0;
        }


        //evaluating Legrange interpolants for each function at xi_point
        lagrange_1D(val_1d, DVal_1d, xi_point(1), nodes_1d, orderN);
        val_2d(m, 1)  = val_1d(j-1); 
        DVal_2d(m, 1) = DVal_1d(j-1);

        // resetting to zero
        for(int i = 0.0; i < N; i++){
            val_1d(i)  = 0.0;
            DVal_1d(i) = 0.0;
        }
        // Assigning and storing the Basis
        lag_basis_2d(m) = val_2d(m, 0) * val_2d(m, 1);

        // Assigning and storing the partials
        lag_partial(m, 0)  = DVal_2d(m, 0) * val_2d(m, 1);
        lag_partial(m, 1)  = val_2d(m, 0) * DVal_2d(m, 1);
    } // end for 

}// end basis_partials function




/* 
 .-------------------------------. 
| .----------------------------. |
| |    ______      ________    | |
| |   / ____ `.   |_   ___ `.  | |
| |   `'  __) |     | |   `. \ | |
| |   _  |__ '.     | |    | | | |
| |  | \____) |    _| |___.' / | |
| |   \______.'   |________.'  | |
| |                            | |
| '----------------------------' |
 '------------------------------' 
*/



/*
==========================
  Hex 8
==========================

 The finite element local vertex numbering for a 8 node Hexahedral is
 as follows

         Mu (k)
         |     Eta (j)    
         |    /
         |   /
     6---+----7
    /|   |   /|
   / |   |  / |
  4--------5  |
  |  |    -|--+---> Xi (i)
  |  |     |  |
  |  2-----|--3
  | /      | /       
  |/       |/
  0--------1
 
*/

real_t Hex8::ref_vert[Hex8::num_verts*Hex8::num_dim] = 
    {// listed as {Xi, Eta, Mu}
    // Bottom Nodes
    -1.0, -1.0, -1.0,// 0
    +1.0, -1.0, -1.0,// 1
    -1.0, +1.0, -1.0,// 2
    +1.0, +1.0, -1.0,// 3
    // Top Nodes
    -1.0, -1.0, +1.0,// 4
    +1.0, -1.0, +1.0,// 5
    -1.0, +1.0, +1.0,// 6
    +1.0, +1.0, +1.0 // 7
    };

const int Hex8::vert_to_node[Hex8::num_verts] = 
    {
    0,
    2,
    6,
    8,
    18,
    20,
    24,
    24
    };

// get the physical location for a given xi_point
void Hex8::physical_position (
    view_c_array <real_t>  &x_point, 
    const view_c_array <real_t>  &xi_point, 
    const view_c_array <real_t>  &vertices){

    real_t basis_a[num_verts];
    auto basis = view_c_array <real_t> (basis_a, num_verts);
    
    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the shape functions from each vertex for (xi,eta,mu)
    for (int vert_lid = 0; vert_lid < num_verts; vert_lid++ ){
        basis(vert_lid) = 1.0/8.0
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // calculate the position in physical space
    for (int dim = 0; dim < num_dim; dim++){
        x_point(dim) = 0.0;
    }

    for (int vert_lid = 0; vert_lid < num_verts; vert_lid++ ){
        for (int dim = 0; dim < num_dim; dim++){
            x_point(dim) += vertices(vert_lid, dim)*basis(vert_lid);
        } // end for dim
    } // end for vert_lid

} // end of function


// calculate the value for the basis at each node for a given xi,eta, mu
void Hex8::basis(
    view_c_array <real_t>  &basis,
    const view_c_array <real_t>  &xi_point){

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the shape functions from each vertex for (xi,eta,mu)
    for (int vert_lid = 0; vert_lid < num_verts; vert_lid++ ){
        basis(vert_lid) = 1.0/8.0
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

} // end of hex8 basis functions


// calculate the partials of the shape function 
// with respect to Xi
void Hex8::partial_xi_basis(
    view_c_array <real_t>  &partial_xi, 
    const view_c_array <real_t>  &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    for (int vert_lid = 0; vert_lid < num_verts; vert_lid++){
        partial_xi(vert_lid) = (1.0/8.0)
            * (ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

} // end of partial Xi function


// with respect to eta
void Hex8::partial_eta_basis(
    view_c_array <real_t> &partial_eta, 
    const view_c_array <real_t>  &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    for (int vert_lid = 0; vert_lid < num_verts; vert_lid++){
        partial_eta(vert_lid) = (1.0/8.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

} //end of partial eta function 


// with repsect to mu
void Hex8::partial_mu_basis(
    view_c_array <real_t> &partial_mu, 
    const view_c_array <real_t> &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    for (int vert_lid = 0; vert_lid < num_verts; vert_lid++){
        partial_mu(vert_lid) = (1.0/8.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (ref_verts(vert_lid, 2));
    } // end for vert_lid

} // end of partial mu function

// Map from vertex to node
inline int Hex8::vert_node_map( const int vert_lid){
    
    return vert_to_node[vert_lid];

};


/*
==========================
  Hex 20
==========================

The finite element local point numbering for a 20 node Hexahedral is 
as follows

           Mu (k)
           |     Eta (j)
           |    /
           |   /

      7----14----6
     /|         /|
   15 |       13 |
   / 19       /  18
  4----12----5   |
  |   |      |   |  --> Xi (i)
  |   |      |   |
  |   3---10-|---2
 16  /      17  /
  | 11       | 9         
  |/         |/
  0-----8----1

*/

real_t Hex20::ref_vert[Hex20::num_verts*Hex20::num_dim] = // listed as {Xi, Eta, Mu}
    // new indices for right hand coordinates
    {
    // bottom corners
    -1.0, -1.0, -1.0,// 0
    +1.0, -1.0, -1.0,// 1
    +1.0, +1.0, -1.0,// 2
    -1.0, +1.0, -1.0,// 3
    // top corners
    -1.0, -1.0, +1.0,// 4
    +1.0, -1.0, +1.0,// 5
    +1.0, +1.0, +1.0,// 6
    -1.0, +1.0, +1.0,// 7
    // bottom edges
     0.0, -1.0, -1.0,// 8
    +1.0,  0.0, -1.0,// 9
     0.0, +1.0, -1.0,// 10
    -1.0,  0.0, -1.0,// 11
    // top edges
     0.0, -1.0, +1.0,// 12
    +1.0,  0.0, +1.0,// 13
     0.0, +1.0, +1.0,// 14
    -1.0,  0.0, +1.0,// 15
    // middle edges
    -1.0, -1.0,  0.0,// 16
    +1.0, -1.0,  0.0,// 17
    +1.0, +1.0,  0.0,// 18
    -1.0, +1.0,  0.0// 19
    };

const int Hex20::vert_to_node[Hex20::num_verts] = 
    {
    0,
    4,
    24,
    20,
    100,
    104,
    124,
    120,
    2,
    14,
    22,
    10,
    102,
    114,
    122,
    110,
    50,
    54,
    74,
    70
    };


// get the physical location for a given xi_point
void Hex20::physical_position (
    view_c_array <real_t>  &x_point, 
    const view_c_array <real_t>  &xi_point, 
    const view_c_array <real_t>  &vertices){

    real_t basis_a[num_verts];
    auto basis = view_c_array <real_t> (basis_a, num_verts);
    
    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the 8 corner shape functions for (xi,eta,mu)
    for (int vert_lid=0; vert_lid<8; vert_lid++){
        basis(vert_lid) = 1.0/8.0
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2))
            * (xi_point(0)*ref_verts(vert_lid, 0)
            +  xi_point(1)*ref_verts(vert_lid, 1)
            +  xi_point(2)*ref_verts(vert_lid, 2) - 2.0);
    } // end for vert_lid

    // calculate the i=0 edge shape functions pts=[8,10,12,14]
    for (int vert_lid = 8; vert_lid <= 14; vert_lid += 2){
        basis(vert_lid) = 1.0/4.0
            * (1.0 - xi_point(0)*xi_point(0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // calculate the k=0 edge shape functions pts=[16,17,18,19]
    for (int vert_lid = 16; vert_lid <= 19; vert_lid++){
        basis(vert_lid) = 1.0/4.0
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 - xi_point(2)*xi_point(2));

   } // end for vert_lid

   // calculate the j=0 edge shape functions pts=[9,11,15,13]
   for (int vert_lid = 9; vert_lid <= 15; vert_lid += 2){
        basis(vert_lid) = 1.0/4.0
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 - xi_point(1)*xi_point(1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // calculate the position in physical space
    for (int dim = 0; dim < num_dim; dim++) x_point(dim) = 0.0;

    for (int dim = 0; dim < num_dim; dim++){
        for (int vert_lid = 0; vert_lid < num_verts; vert_lid++ ){
            x_point(dim) += vertices(vert_lid, dim)*basis(vert_lid);
        }   
    } // end for dim

} // end of physical position function


// calculate the value for the basis at each node for a given xi,eta, mu
void Hex20::basis(
    view_c_array <real_t>  &basis,
    const view_c_array <real_t>  &xi_point){

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the 8 corner shape functions for (xi,eta,mu)
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
        basis(vert_lid) = 1.0/8.0
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2))
            * (xi_point(0)*ref_verts(vert_lid, 0)
            +  xi_point(1)*ref_verts(vert_lid, 1)
            +  xi_point(2)*ref_verts(vert_lid, 2) - 2.0);
    } // end for vert_lid

    // calculate the i=0 edge shape functions pts=[8,10,12,14]
    for (int vert_lid = 8; vert_lid <= 14; vert_lid += 2){
        basis(vert_lid) = 1.0/4.0
            * (1.0 - xi_point(0)*xi_point(0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // calculate the k=0 edge shape functions pts=[16,17,18,19]
    for (int vert_lid = 16; vert_lid <= 19; vert_lid++){
        basis(vert_lid) = 1.0/4.0
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 - xi_point(2)*xi_point(2)); 
    } // end for vert_lid

    // calculate the j=0 edge shape functions pts=[9,11,15,13]
    for (int vert_lid = 9; vert_lid <= 15; vert_lid += 2){
        basis(vert_lid) = 1.0/4.0
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 - xi_point(1)*xi_point(1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));           
    } // end for vert_lid

} // end of hex20 basis functions


// Calculate the partials of the shape functions
// with respect to Xi
void  Hex20::partial_xi_basis(
    view_c_array <real_t>  &partial_xi, 
    const view_c_array <real_t>  &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // For 8 Corner shape functions pts=[0,1,2,3,4,5,6,7]
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
        partial_xi(vert_lid) = (1.0/8.0) 
            * (ref_verts(vert_lid, 0))
            * (1.0 + (xi_point(1)*ref_verts(vert_lid, 1)))
            * (1.0 + (xi_point(2)*ref_verts(vert_lid, 2)))
            * (2.0 * (xi_point(0)*ref_verts(vert_lid, 0))
            + xi_point(1)*ref_verts(vert_lid, 1)
            + xi_point(2)*ref_verts(vert_lid, 2) - 1.0);  
    } // end for vert_lid

    // for the i=0 edge shape functions pts=[8,10,12,14]
    for (int vert_lid = 8; vert_lid <= 14; vert_lid += 2){
        partial_xi(vert_lid) = (-1.0/2.0)
            * (xi_point(0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // for the k=0 shape functions pts=[9,11,13,15]
    for (int vert_lid = 9; vert_lid <= 15; vert_lid += 2){
        partial_xi(vert_lid) = (1.0/4.0)
            * (ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 - xi_point(2)*xi_point(2));
    } // end for vert_lid


    // for the j=0 edge shape functions pts=[16,17,18,19]
    for (int vert_lid = 16; vert_lid <= 19; vert_lid++){
        partial_xi(vert_lid) = (1.0/4.0)
            * (ref_verts(vert_lid, 0))
            * (1.0 - xi_point(1)*xi_point(1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

} // end of partial Xi function


// with respect to Eta
void Hex20::partial_eta_basis(
    view_c_array <real_t> &partial_eta, 
    const view_c_array <real_t>  &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // For 8 Corner shape functions pts=[0,1,2,3,4,5,6,7]
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
        partial_eta(vert_lid) = (1.0/8.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2))
            * (xi_point(0)*ref_verts(vert_lid, 0)
            +  2.0 * xi_point(1)*ref_verts(vert_lid, 1)
            + xi_point(2)*ref_verts(vert_lid, 2) - 1.0);
    } // end for vert_lid

    // for the i=0 edge shape functions pts=[8,10,12,14]
    for (int vert_lid = 8; vert_lid <= 14; vert_lid += 2){
        partial_eta(vert_lid) = (1.0/4.0)
            * (1.0 - (xi_point(0)*xi_point(0)))
            * (ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // for the j=0 shape functions pts=[9,11,13,15]
    for (int vert_lid = 9; vert_lid <= 15; vert_lid += 2){
        partial_eta(vert_lid) = (-1.0/2.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (xi_point(1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // for the k=0 edge shape functions pts=[16,17,18,19]
    for (int vert_lid = 16; vert_lid <= 19; vert_lid++){
        partial_eta(vert_lid) = (1.0/4.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (ref_verts(vert_lid, 1))
            * (1.0 - (xi_point(2)*xi_point(2)));
    } // end for vert_lid

} // end of partial Eta function


// with repsect to mu
void Hex20::partial_mu_basis(
    view_c_array <real_t> &partial_mu, 
    const view_c_array <real_t>  &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // For 8 Corner shape functions pts=[0,1,2,3,4,5,6,7]
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
        partial_mu(vert_lid) = (1.0/8.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (ref_verts(vert_lid, 2))
            * ((xi_point(0)*ref_verts(vert_lid, 0))
            + (xi_point(1)*ref_verts(vert_lid, 1))
            + (2.0 * xi_point(2)*ref_verts(vert_lid, 2)) - 1.0);
    } // end for vert_lid

    // for the i=0 edge shape functions pts=[8,10,12,14]
    for (int vert_lid = 8; vert_lid <= 14; vert_lid += 2){
        partial_mu(vert_lid) = (1.0/4.0)
            * (1.0 - (xi_point(0)*xi_point(0))) 
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (ref_verts(vert_lid, 2));
    }

    // for the j=0 shape functions pts=[9,11,13,15]
    for (int vert_lid = 9; vert_lid <= 15; vert_lid += 2){
        partial_mu(vert_lid) = (1.0/4.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 - xi_point(1)*xi_point(1))
            * (ref_verts(vert_lid, 2));
    } // end for vert_lid

    // for the j=0 edge shape functions pts=[16,17,18,19]
    for (int vert_lid = 16; vert_lid <= 19; vert_lid++){
        partial_mu(vert_lid) = (-1.0/2.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (xi_point(2));
    } // end for vert_lid

} // end of partial Mu function



// Map from vertex to node
inline int Hex20::vert_node_map( const int vert_lid){

    return vert_to_node[vert_lid];
};

/* 
==========================
  Hex 32
==========================

The finite element local point numbering for a 32 node Hexahedral is 
shown below


               Mu (k)
                ^         Eta (j)
                |        /
                |       /
                       /
        7----23------22----6
       /|                 /|
     15 |               14 |
     /  |               /  |
   12  31             13   30 
   /    |             /    |
  4-----20-----21----5     |
  |     |            |     |   ----> Xi (i)
  |    27            |     26  
  |     |            |     |
 28     |           29     |
  |     3----19------|18---2
  |    /             |    /
  |  11              |   10
 24  /              25  /
  | 8                | 9         
  |/                 |/
  0----16------17----1
*/




real_t Hex32::ref_vert[Hex32::num_verts*Hex32::num_dim] = // listed as {Xi, Eta, Mu}
    {
    -1.0, -1.0, -1.0,// 0
    +1.0, -1.0, -1.0,// 1
    +1.0, +1.0, -1.0,// 2
    -1.0, +1.0, -1.0,// 3
    -1.0, -1.0, +1.0,// 4
    +1.0, -1.0, +1.0,// 5
    +1.0, +1.0, +1.0,// 6
    -1.0, +1.0, +1.0,// 7
    // Xi/Mu = +- 1/3
    -1.0, -1./3., -1.0,// 8
     1.0, -1./3., -1.0,// 9
     1.0, +1./3., -1.0,// 10
    -1.0, +1./3., -1.0,// 11
    -1.0, -1./3., +1.0,// 12
     1.0, -1./3., +1.0,// 13
     1.0, +1./3., +1.0,// 14
    -1.0, +1./3., +1.0,// 15
    // Eta/Mu = +- 1/3
    -1./3., -1.0, -1.0,// 16
     1./3., -1.0, -1.0,// 17
     1./3., +1.0, -1.0,// 18
    -1./3., +1.0, -1.0,// 19
    -1./3., -1.0,  1.0,// 20
     1./3., -1.0,  1.0,// 21
     1./3., +1.0,  1.0,// 22
    -1./3., +1.0,  1.0,// 23
    // Xi/Eta = +- 1/3
    -1.0, -1.0, -1./3.,// 24
     1.0, -1.0, -1./3.,// 25
     1.0,  1.0, -1./3.,// 26
    -1.0,  1.0, -1./3.,// 27
    -1.0, -1.0,  1./3.,// 28
     1.0, -1.0,  1./3.,// 29
     1.0,  1.0,  1./3.,// 30
    -1.0,  1.0,  1./3.,// 31
    };

const int Hex32::vert_to_node[Hex32::num_verts] = 
    {
    0,
    6,
    48,
    42,
    294,
    300,
    342,
    336,
    14,
    20,
    32,
    28,
    308,
    314,
    328,
    322,
    2,
    4,
    46,
    44,
    296,
    298,
    340,
    338,
    98,
    104,
    146,
    140,
    196,
    202,
    244,
    298
    };

// get the physical location for a given xi_point
void Hex32::physical_position (
    view_c_array <real_t>  &x_point, 
    const view_c_array <real_t>  &xi_point, 
    const view_c_array <real_t>  &vertices){

    real_t basis_a[num_verts];
    auto basis = view_c_array <real_t> (basis_a, num_verts);
    
    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the 8 corner shape functions for (xi,eta,mu)
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
        basis(vert_lid) = (1.0/64.0) 
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2))
            * (9.0 * xi_point(0)*xi_point(0)
            +  9.0 * xi_point(1)*xi_point(1)
            +  9.0 * xi_point(2)*xi_point(2) - 19.0);  
    } // end for vert_lid

    // calculate the edge shape functions for pts=[8-15]
    for (int vert_lid = 8; vert_lid <= 15; vert_lid++){
        basis(vert_lid) = (9.0/64.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 - xi_point(1)*xi_point(1))
            * (1.0 + 9.0 * xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid
    
    // calculate the edge shape functions for pts=[16-23]
    for (int vert_lid = 16; vert_lid <= 23; vert_lid++){
        basis(vert_lid) = (9.0/64.0)
            * (1.0 - xi_point(0)*xi_point(0))
            * (1.0 + 9.0 * xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2)); 
    } // end for vert_lid

    // calculate the edge shape functions for pts=[24-31]
    for (int vert_lid = 24; vert_lid <= 31; vert_lid++){
        basis(vert_lid) = (9.0/64.0) 
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + 9.0 * xi_point(2)*ref_verts(vert_lid, 2))
            * (1.0 - xi_point(2)*xi_point(2)); 
    } // end for vert_lid


    // calculate the position in physical space
    for (int dim = 0; dim < num_dim; dim++){
        x_point(dim) = 0.0;
    }

    for (int vert_lid = 0; vert_lid <= num_verts; vert_lid++ ){
        //std::cout << "Vert :" << vert_lid << std::endl;
        for (int dim = 0; dim < num_dim; dim++){
            x_point(dim) += vertices(vert_lid, dim)*basis(vert_lid);
        }
    }

} // end of physical position function


void Hex32::basis(
    view_c_array <real_t>  &basis,
    const view_c_array <real_t>  &xi_point){

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the 8 corner shape functions for (xi,eta,mu)
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
        basis(vert_lid) = (1.0/64.0) 
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2))
            * (9.0 * xi_point(0)*xi_point(0)
            +  9.0 * xi_point(1)*xi_point(1)
            +  9.0 * xi_point(2)*xi_point(2) - 19.0);  
    } // end for vert_lid

    // calculate the edge shape functions for pts=[8-15]
    for (int vert_lid = 8; vert_lid <= 15; vert_lid++){
        basis(vert_lid) = (9.0/64.0)
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 - xi_point(1)*xi_point(1))
            * (1.0 + 9.0 * xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // calculate the edge shape functions for pts=[16-23]
    for (int vert_lid = 16; vert_lid <= 23; vert_lid++){
        basis(vert_lid) = (9.0/64.0)
            * (1.0 - xi_point(0)*xi_point(0))
            * (1.0 + 9.0 * xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2)*ref_verts(vert_lid, 2)); 
    } // end for vert_lid

    // calculate the edge shape functions for pts=[24-31]
    for (int vert_lid = 24; vert_lid < num_verts; vert_lid++){
        basis(vert_lid) = (9.0/64.0) 
            * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
            * (1.0 + 9.0 * xi_point(2)*ref_verts(vert_lid, 2))
            * (1.0 - xi_point(2)*xi_point(2)); 
    } // end for vert_lid

} // end of hex20 basis functions

// Calculate the partials of the shape functions
// with respect to Xi
void  Hex32::partial_xi_basis(
    view_c_array <real_t>  &partial_xi, 
    const view_c_array <real_t>  &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the 8 corner partial wrt Xi 
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
        partial_xi(vert_lid) = (1.0/64.0) 
            * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2) * ref_verts(vert_lid, 2))
            *((9.0 * (ref_verts(vert_lid, 0))
            * (xi_point(0)*xi_point(0) +  xi_point(1)*xi_point(1) + xi_point(2)*xi_point(2)))
            + (18.0 * xi_point(0) * (1.0 + xi_point(0)*ref_verts(vert_lid, 0)))
            - (19.0 * ref_verts(vert_lid, 0)));
    } // end for vert_lid

    // calculate the edge partial wrt Xi for pts=[8-15]
    for (int vert_lid = 8; vert_lid <= 15; vert_lid++){
        partial_xi(vert_lid) = (9.0/64.0) 
            * (ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
            * (1.0 + 9.0 * xi_point(2) * ref_verts(vert_lid, 2))
            * (1.0 - xi_point(2) * xi_point(2));
    } // end for vert_lid

    // calculate the edge partial wrt Xi for pts=[16-23]
    for (int vert_lid = 16; vert_lid <= 23; vert_lid++){
        partial_xi(vert_lid) = (9.0/64.0) 
            * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2) * ref_verts(vert_lid, 2))
            * (9.0 * ref_verts(vert_lid, 0) * (1.0 - 3.0 * xi_point(0) * xi_point(0))
            - (2 * xi_point(0)));
    } // end for vert_lid

    // calculate the edge partial wrt Xi for pts=[24-31]
    for (int vert_lid = 24; vert_lid <= 31; vert_lid++){
        partial_xi(vert_lid) = (9.0/64.0) 
            * (ref_verts(vert_lid, 0))
            * (1.0 - xi_point(1) * xi_point(1))
            * (1.0 + 9.0 * xi_point(1) * ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2) * ref_verts(vert_lid, 2));
    } // end for vert_lid

} // end of partial Xi function


// with respect to Eta
// functions for [18-15] and [24-31] were switched 
void Hex32::partial_eta_basis(
    view_c_array <real_t> &partial_eta, 
    const view_c_array <real_t>  &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the 8 corner partial wrt Eta 
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
        partial_eta(vert_lid) = (1.0/64.0) 
            * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
            * (1.0 + xi_point(2) * ref_verts(vert_lid, 2))
            *((9.0 * ref_verts(vert_lid, 1)
            * (xi_point(0)*xi_point(0) 
            +  xi_point(1)*xi_point(1)
            +  xi_point(2)*xi_point(2)))
            + (18.0 * xi_point(1) * (1.0 + xi_point(1)*ref_verts(vert_lid, 1)))
            - (19.0 * ref_verts(vert_lid, 1)));
    } // end for vert_lid

    // calculate the edge partial wrt Eta for pts=[8-15]
    for (int vert_lid = 8; vert_lid <= 15; vert_lid++){
        partial_eta(vert_lid) = (9.0/64.0) 
            * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
            * (1.0 + xi_point(2) * ref_verts(vert_lid, 2))
            *((9.0 * ref_verts(vert_lid, 1) * (1.0 - 3.0 * xi_point(1) * xi_point(1)))
            - (2.0 * xi_point(1)));
    } // end for vert_lid

    // calculate the edge partial wrt Eta for pts=[16-23]
    for (int vert_lid = 16; vert_lid <= 23; vert_lid++){
        partial_eta(vert_lid) = (9.0/64.0) 
            * (1.0 - xi_point(0) * xi_point(0))
            * (1.0 + 9.0 * xi_point(0) * ref_verts(vert_lid, 0))
            * (ref_verts(vert_lid, 1))
            * (1.0 + xi_point(2) * ref_verts(vert_lid, 2));
    } // end for vert_lid

    // calculate the edge partial wrt Eta for pts=[24-31]
    for (int vert_lid = 24; vert_lid <= 31; vert_lid++){
        partial_eta(vert_lid) = (9.0/64.0) 
            * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
            * (ref_verts(vert_lid, 1))
            * (1.0 + 9.0 * xi_point(2) * ref_verts(vert_lid, 2))
            * (1.0 - xi_point(2) * xi_point(2));
    } // end for vert_lid

} // end of partial Eta function


// with repsect to mu
// functions for [18-15] and [24-31] were switched 
void Hex32::partial_mu_basis(
    view_c_array <real_t> &partial_mu, 
    const view_c_array <real_t>  &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the 8 corner partial wrt Mu 
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
        partial_mu(vert_lid) = (1.0/64.0)
            * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
            *((9.0 * (ref_verts(vert_lid, 2))
            * (xi_point(0)*xi_point(0) 
            +  xi_point(1)*xi_point(1)
            +  xi_point(2)*xi_point(2))) 
            + (18.0 * xi_point(2) * (1 + xi_point(2)*ref_verts(vert_lid, 2)))
            - (19.0 * ref_verts(vert_lid, 2)));
    } // end for vert_lid

    // calculate the edge partial wrt Mu for pts=[8-15]
    for (int vert_lid = 8; vert_lid <= 15; vert_lid++){
        partial_mu(vert_lid) = (9.0/64.0) 
            * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
            * (1.0 - xi_point(1) * xi_point(1))
            * (1.0 + 9.0 * xi_point(1) * ref_verts(vert_lid, 1))
            * (ref_verts(vert_lid, 2));

    } // end for vert_lid

    // calculate the edge partial wrt Mu for pts=[16-23]
    for (int vert_lid = 16; vert_lid <= 23; vert_lid++){
        partial_mu(vert_lid) = (9.0/64.0) 
            * (1.0 - xi_point(0) * xi_point(0))
            * (1.0 + 9.0 * xi_point(0) * ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
            * (ref_verts(vert_lid, 2));
    } // end for vert_lid

    // calculate the edge partial wrt Mu for pts=[24-31]
    for (int vert_lid = 24; vert_lid <= 31; vert_lid++){
        partial_mu(vert_lid) = (9.0/64.0) 
            * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
            * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
            *((9.0 * ref_verts(vert_lid, 2) 
            * (1.0 -  3.0 * xi_point(2) * xi_point(2)))
            - (2.0 * xi_point(2)));
    } // end for vert_lid

} // end of partial Mu function


/* 
==========================
  Hex 32
==========================

The finite element local point numbering for a 32 node Hexahedral is 
shown below


               Mu (k)
                ^         Eta (j)
                |        /
                |       /
                       /
        7----23------22----6
       /|                 /|
     15 |               14 |
     /  |               /  |
   12  31             13   30 
   /    |             /    |
  4-----20-----21----5     |
  |     |            |     |   ----> Xi (i)
  |    27            |     26  
  |     |            |     |
 28     |           29     |
  |     3----19------|18---2
  |    /             |    /
  |  11              |   10
 24  /              25  /
  | 8                | 9         
  |/                 |/
  0----16------17----1
*/

// Map from vertex to node
inline int Hex32::vert_node_map( const int vert_lid){

    return vert_to_node[vert_lid];
};





/*
  _   _           _   _ 
 | | | | _____  _| \ | |
 | |_| |/ _ \ \/ /  \| |
 |  _  |  __/>  <| |\  |
 |_| |_|\___/_/\_\_| \_|
                        
representative linear element for visualization
   
            j
            |     k    
            |    /
            |   /
        6---+----7
       /|   |   /|
      / |   |  / |
     2--------3  |
     |  |    -|--+---> i
     |  |     |  |
     |  4-----|--5
     | /      | /       
     |/       |/
     0--------1
    
   Note: left hand coordinate coordinates

*/


// Lagrange Interp in 1D, returns interpolants and derivative
// works with any nodal spacing
void HexN::lagrange_1D(
    view_c_array <real_t> &interp,          // interpolant
    view_c_array <real_t> &Dinterp,         // derivative of function
    const real_t &x_point,                  // point of interest in element
    const view_c_array <real_t> &xi_point,  // nodal positions in 1D, normally chebyshev
    const int &orderN){                     // order of element

    real_t num_a[orderN+1];
    auto num = view_c_array <real_t> (num_a, orderN+1); // numerator of interpolant

    real_t denom_a[orderN+1];
    auto denom = view_c_array <real_t> (denom_a, orderN+1); // denomenator of interpolant
  
    real_t q = 0.0;
   
    for(int i = 0; i < orderN + 1; i++){ // looping over the nodes
        real_t n = 1.0;         // placeholder numerator
        real_t d = 1.0;         // placeholder denominator
        real_t c = 1.0;         // placeholder value of n/d
        real_t p = 0.0;         // placeholder for derivative values
        real_t s = 1.0;
        
        for(int j = 0; j < orderN + 1; j++){  // looping over the nodes !=i
            if (j != i ){
                n = n*(x_point - xi_point(j));
                d = d*(xi_point(i) - xi_point(j));
                real_t s = 1.0;
                
                for(int N = 0; N < orderN + 1; N++){  // looping over the nodes !=i
                    if (N != j && N != i ){
                    s = s * (x_point - xi_point(N));
                    }// end if
                }//end for
                
                p += s; 
            }//end if
         
            c = n/d; // storing a single value for interpolation for node i
            q = (p/d); // storing the derivative of the interpolating function
        } // end looping over nodes != i

        // writing value to vectors for later use
        interp(i)   = c;     // Interpolant value at given point
        Dinterp(i)  = q;     // derivative of each function
    } // end loop over all nodes

} // end of Legrange_1D function


// Corners of Lagrange element for mapping
void HexN::corners (
    view_c_array <real_t> &lag_nodes,   // Nodes of Lagrange elements 
    view_c_array <real_t> &lag_corner,  // corner nodes of HexN element
    const int &orderN){                 // Element order


   /*
   This image represents the corner mapping notation of an arbitrary ordered
   Lagrange element. The corner function takes in the element order and nodal positions and
   returns a vector containing the indices of the corner in alphabetical order.

           j
           |     k    
           |    /
           |   /
       G---+----H
      /|   |   /|
     / |   |  / |
    C--------D  |
    |  |    -|--+---> i
    |  |     |  |
    |  E-----|--F
    | /      | /       
    |/       |/
    A--------B
   
   Note: left hand coordinate coordinates
   */
    int num_corners = 8;
    int N = orderN + 1;      //number of nodes in each direction
    int corner_ids[num_corners];
    

    corner_ids[0] = 0;                      
    corner_ids[1] = N - 1.0;                 
    corner_ids[2] = (N*N) - N;         
    corner_ids[3] = (N*N)-1.0;               
    corner_ids[4] = (N*N*N) - (N*N) ;  
    corner_ids[5] = (N*N*N) - (N*N) + N - 1.0;  
    corner_ids[6] = (N*N*N) - (N);     
    corner_ids[7] = (N*N*N) - 1; 

    for(int corner = 0; corner < num_corners; corner++){
        for(int dim = 0; dim < num_dim; dim++){
            lag_corner(corner, dim) = lag_nodes(corner_ids[corner], dim);
        }
    }

    // lag_corner[0] = lag_nodes[A];
    // std::cout<<"Corner A = "<<lag_corner[0][0] << "  "<<lag_corner[0][1] << "  "<<lag_corner[0][2] <<std::endl;

    // lag_corner[1] = lag_nodes[B];
    // std::cout<<"Corner B = "<<lag_corner[1][0] << "  "<<lag_corner[1][1] << "  "<<lag_corner[1][2] <<std::endl;

    // lag_corner[2] = lag_nodes[C];
    // std::cout<<"Corner C = "<<lag_corner[2][0] << "  "<<lag_corner[2][1] << "  "<<lag_corner[2][2] <<std::endl;

    // lag_corner[3] = lag_nodes[D];
    // std::cout<<"Corner D = "<<lag_corner[3][0] << "  "<<lag_corner[3][1] << "  "<<lag_corner[3][2] <<std::endl;

    // lag_corner[4] = lag_nodes[E];
    // std::cout<<"Corner E = "<<lag_corner[4][0] << "  "<<lag_corner[4][1] << "  "<<lag_corner[4][2] <<std::endl;

    // lag_corner[5] = lag_nodes[F];
    // std::cout<<"Corner F = "<<lag_corner[5][0] << "  "<<lag_corner[5][1] << "  "<<lag_corner[5][2] <<std::endl;

    // lag_corner[6] = lag_nodes[G];
    // std::cout<<"Corner G = "<<lag_corner[6][0] << "  "<<lag_corner[6][1] << "  "<<lag_corner[6][2] <<std::endl;

    // lag_corner[7] = lag_nodes[H];
    // std::cout<<"Corner H = "<<lag_corner[7][0] << "  "<<lag_corner[7][1] << "  "<<lag_corner[7][2] <<std::endl;

}// end of corner mapping function


// Functions for mapping reference position to physical position for any 
// point in an arbitrary order 3D lagrange element
void HexN::physical_position (
    view_c_array <real_t> &x_point,             // location in real space
    const view_c_array <real_t> &lag_nodes,     // Nodes of Lagrange elements 
    const view_c_array <real_t> &lag_basis_3d,  // 3D basis values 
    const int &orderN){                         // order of the element

    int nodes = orderN + 1;
    int Nnodes_3d = nodes * nodes * nodes;

    for (int this_vert = 0; this_vert < Nnodes_3d; this_vert++ ){
        for (int dim = 0; dim < num_dim; dim++){
            x_point(dim) += lag_nodes(this_vert, dim)*lag_basis_3d(this_vert);
        } // end for dim
    } // end for this_vert

}// end physical position function


void HexN::basis_partials (
    view_c_array <real_t> &lag_nodes,       // Nodes of Lagrange elements (to be filled in)
    view_c_array <real_t> &nodes_1d,        // Nodal spacing in 1D, any spacing is accepted
    view_c_array <real_t> &val_1d,          // Interpolant Value in 1D
    view_c_array <real_t> &DVal_1d,         // Derivateive of basis in 1D
    view_c_array <real_t> &val_3d,          // for holding the interpolant in each direction
    view_c_array <real_t> &DVal_3d,         // for holding the derivatives in each direction
    view_c_array <real_t> &lag_basis_3d,    // 3D basis values 
    view_c_array <real_t> &lag_partial,     // Partial of basis 
    const view_c_array <real_t> &xi_point,  // point of interest
    const int &orderN){                     // Element order

   /*

   representative linear element for visualization
   
            j
            |     k    
            |    /
            |   /
        6---+----7
       /|   |   /|
      / |   |  / |
     2--------3  |
     |  |    -|--+---> i
     |  |     |  |
     |  4-----|--5
     | /      | /       
     |/       |/
     0--------1
    
   Note: left hand coordinate coordinates
   */

    int N = orderN + 1;      //number of nodes in each direction
    int tot_pts = (N*N*N);  // total nodes in 3D

    real_t sumi = 0.0;
    real_t sumj = 0.0;
    real_t sumk = 0.0;


    //Setting nodes for Lagrange Elements
    for (int m = 0; m < tot_pts; m++) {

        int i, j, k, jj;
        // sets up the i, j, k indices for the nodes of an 
        // arbitrary order lagrange element
        int p12 = (N * N);

        k  = floor(m/p12) + 1; 

        jj = (m+1) - p12 * (k-1); 
        j  = floor((jj-1) / N)+1;

        i  = jj - N * (j-1); 

        // xi direction
        lag_nodes(m, 0) = nodes_1d(i-1); 

        // eta direction
        lag_nodes(m, 1) = nodes_1d(j-1); 

        // mu direction 
        lag_nodes(m, 2) = nodes_1d(k-1); 

        // calling function to assign nodal values for basis and derivative

        //evaluating Lagrange interpolants for each function at xi_point[0]
        lagrange_1D(val_1d, DVal_1d, xi_point(0), nodes_1d, orderN);
        val_3d(m, 0)  = val_1d(i-1); 
        DVal_3d(m, 0) = DVal_1d(i-1);

        // resetting to zero
        for(int i = 0.0; i < N; i++){
            val_1d(i)  = 0.0;
            DVal_1d(i) = 0.0;
        }


        //evaluating Legrange interpolants for each function at xi_point[1]
        lagrange_1D(val_1d, DVal_1d, xi_point(1), nodes_1d, orderN);
        val_3d(m, 1)  = val_1d(j-1); 
        DVal_3d(m, 1) = DVal_1d(j-1);

        // resetting to zero
        for(int i = 0.0; i < N; i++){
            val_1d(i)  = 0.0;
            DVal_1d(i) = 0.0;
        }

        //evaluating Legrange interpolants for each function at xi_point[2]
        lagrange_1D(val_1d, DVal_1d, xi_point(2), nodes_1d, orderN);
        val_3d(m, 2)  = val_1d(k-1);   //3d[2dto1d(m, 2)]
        DVal_3d(m, 2) = DVal_1d(k-1);
        // resetting to zero
        for(int i = 0.0; i < N; i++){
            val_1d(i)  = 0.0;
            DVal_1d(i) = 0.0;
        }


        // Assigning and storing the Basis
        lag_basis_3d(m) = val_3d(m, 0) * val_3d(m, 1) * val_3d(m, 2);

        // Assigning and storing the partials

        lag_partial(m, 0)  = DVal_3d(m, 0) * val_3d(m, 1) * val_3d(m, 2);
        lag_partial(m, 1)  = val_3d(m, 0) * DVal_3d(m, 1) * val_3d(m, 2);
        lag_partial(m, 2)  = val_3d(m, 0) * val_3d(m, 1) * DVal_3d(m, 2);
    } // end for  

}// end basis_partials function




/*
==========================
   4D Tesseract element
==========================

The finite element local point numbering for a 16 node Tesseract is
based on the 3D Hex8 Ensight element
 

                 _.15-------------------------------------14
            _.+<    |\                              . >-"/ |
      _ .+>         | \                         .>"" ./    |
  .>""              |  \                     <""    /      |
12----------------------+------------------13    ./        |
| )<=               |    \               / | _/""          |
|     )\+           |     \            / /"|               |
|         (\=       |   _. 7---------+--6  |               |
|             \>   .|+<    |       / . "|  |               |
|               '4--+------+------5'    |  |               |
|                |  |      |      |     |  |               |
|                |  |      |      |     |  |               |
|                |  |      |      |     |  |               |
|                |  |      |      |     |  |               |
|                |  |      |      |     |  |               |
|                |  |      |      |     |  |               |
|                |  |   _ .3------+-----2_ |               |
|                |  | "   /       |   /'  '| \= _          |
|                0--+---+---------1*"      |     ""\       |
|             ./'   |  "           \       |         "">   |
|            /      |/              \      |             ".|
|         /        11----------------+-----+--------------10
|      ./    .+<""                    )    |           .</
|    /(   /(                            \  |     _.+</
| ./  /"                                 \ |  >(
8------------------------------------------9'

              j
              ^        k
              |      /
              |    / 
              |  /
              |/
              +---------->i

            i = Xi
            j = Eta
            k = Mu
            t = Tau 

*/




real_t Tess16::ref_vert[ Tess16::num_verts* Tess16::num_dim] = // listed as {Xi, Eta, Mu, Tau}
    {
    // Interior cube bottom
    -1.0, -1.0, -1.0, -1.0,
    +1.0, -1.0, -1.0, -1.0,
    +1.0, -1.0, +1.0, -1.0,
    -1.0, -1.0, +1.0, -1.0,
    // Interior cube top
    -1.0, +1.0, -1.0, -1.0,
    +1.0, +1.0, -1.0, -1.0,
    +1.0, +1.0, +1.0, -1.0,
    -1.0, +1.0, +1.0, -1.0,
    // Exterior cube bottom
    -1.0, -1.0, -1.0, +1.0,
    +1.0, -1.0, -1.0, +1.0,
    +1.0, -1.0, +1.0, +1.0,
    -1.0, -1.0, +1.0, +1.0,
    // Exterior cube top
    -1.0, +1.0, -1.0, +1.0,
    +1.0, +1.0, -1.0, +1.0,
    +1.0, +1.0, +1.0, +1.0,
    -1.0, +1.0, +1.0, +1.0,
    };


// calculate a physical position in an element for a given xi,eta,mu
void Tess16::physical_position(
    view_c_array <real_t> &x_point,
    const view_c_array <real_t> &xi_point,
    const view_c_array <real_t> &vertices){

    real_t basis_a[num_verts];
    auto basis = view_c_array <real_t> (basis_a, num_verts);
    
    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);
   
    // calculate the shape functions from each vertex for (xi,eta,mu, tau)
    for(int this_vert = 0; this_vert < num_verts; this_vert++){
        basis(this_vert) = 1.0/16.0
            * (1.0 + xi_point(0)*ref_verts(this_vert, 0)) 
            * (1.0 + xi_point(1)*ref_verts(this_vert, 1)) 
            * (1.0 + xi_point(2)*ref_verts(this_vert, 2)) 
            * (1.0 + xi_point(3)*ref_verts(this_vert, 3));
    } // end for shape functions

    // calculate the position in physical space
    for (int dim = 0; dim < 4; dim++){
        x_point(dim) = 0.0;
    }

    for (int this_vert = 0; this_vert < 16; this_vert++ ){
        for (int dim = 0; dim < 4; dim++){
            x_point(dim) += vertices(this_vert, dim)*basis(this_vert);
        } // end for dim
    } // end for this_vert

} // End physical position function


// calculate the value for the basis at each node for a given xi,eta,mu,tau
void Tess16::basis(
    view_c_array <real_t>  &basis,
    const view_c_array <real_t>  &xi_point){

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    // calculate the basis functions from each vertex for (xi,eta,mu, tau)
    for(int this_vert = 0; this_vert < num_verts; this_vert++){
        basis(this_vert) = 1.0/16.0
            * (1.0 + xi_point(0)*ref_verts(this_vert, 0)) 
            * (1.0 + xi_point(1)*ref_verts(this_vert, 1)) 
            * (1.0 + xi_point(2)*ref_verts(this_vert, 2)) 
            * (1.0 + xi_point(3)*ref_verts(this_vert, 3));
    } // end for this_vert

}


// Partial derivative of shape functions with respect to Xi at Xi_point
void Tess16::partial_xi_basis(
    view_c_array <real_t> &partial_xi, 
    const view_c_array <real_t> &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    for (int this_vert = 0; this_vert < num_verts; this_vert++){
        partial_xi(this_vert) = 1.0/16.0
            * (ref_verts(this_vert, 0))
            * (1.0 + xi_point(1)*ref_verts(this_vert, 1))
            * (1.0 + xi_point(2)*ref_verts(this_vert, 2))
            * (1.0 + xi_point(3)*ref_verts(this_vert, 3));
    }   // end for this_vert 

} // end partial Xi function


// Partial derivative of shape functions with respect to Eta
void Tess16::partial_eta_basis(
    view_c_array <real_t> &partial_eta, 
    const view_c_array <real_t> &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    for (int this_vert = 0; this_vert < num_verts; this_vert++){  
        partial_eta(this_vert) = 1.0/16.0
            * (1.0 + xi_point(0)*ref_verts(this_vert, 0))
            * (ref_verts(this_vert, 1))
            * (1.0 + xi_point(2)*ref_verts(this_vert, 2))
            * (1.0 + xi_point(3)*ref_verts(this_vert, 3));               
    }   // end for this_vert 

}  // End partial eta function


// Partial derivative of shape functions with respect to Mu
void Tess16::partial_mu_basis(
    view_c_array <real_t> &partial_mu, 
    const view_c_array <real_t> &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    for (int this_vert = 0; this_vert < num_verts; this_vert++){  
        partial_mu(this_vert) = 1.0/16.0
            * (1.0 + xi_point(0)*ref_verts(this_vert, 0))
            * (1.0 + xi_point(1)*ref_verts(this_vert, 1))
            * (ref_verts(this_vert, 2))
            * (1.0 + xi_point(3)*ref_verts(this_vert, 3));
    } // end for this_vert 

} // end partial Mu fuction


// Partial derivative of shape functions with respect to Tau
void Tess16::partial_tau_basis(
    view_c_array <real_t> &partial_tau, 
    const view_c_array <real_t> &xi_point) {

    auto ref_verts = view_c_array<real_t> (ref_vert, num_verts, num_dim);

    for (int this_vert = 0; this_vert < num_verts; this_vert++){  
        partial_tau(this_vert) = 1.0/16.0
            * (1.0 + xi_point(0)*ref_verts(this_vert, 0))
            * (1.0 + xi_point(1)*ref_verts(this_vert, 1))
            * (1.0 + xi_point(2)*ref_verts(this_vert, 2))
            * (ref_verts(this_vert, 3));
    } // end for this_vert   

} // End partial tau function


} // end namespace elements::swage
} // end namespace elements
