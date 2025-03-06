Shader "Unlit/CrRGB"
{
    Properties
    {
        _ToneImg ("Texture", 2D) = "white" {}

        _InvGamma ("Inv Gamma", Color) = (0.,0.,0.)
        _BlRGB ("Bl RGB", Color) = (0.,0.,0.)
        _MaxLum ("Max Lum", Float) = 0.0
        _FlickerCount ("Flicker Count", int) = 0
        _GamutWarning ("Gamut Warning", int) = 0
        _HWBitDepth ("HW Bit Depth", int) = 0

    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

      Pass { // some shaders require multiple passes
         GLSLPROGRAM // here begins the part in Unity's GLSL
        #pragma multi_compile_local _ ADJUST_TO_LINEARSPACE


        #extension GL_ARB_texture_rectangle : enable
        #extension GL_EXT_gpu_shader4 : enable

        uniform sampler2D _ToneImg;
        uniform vec3 _InvGamma;


        uniform mat3x3 _M_rec709_rgbdisp;
        uniform vec3 _BlRGB;
        uniform float _MaxLum;

        uniform highp int _FlickerCount;
        uniform highp int _GamutWarning;
        uniform float _HWBitDepth;


         #ifdef VERTEX // here begins the vertex shader

         out vec4 texCoord;
         void main() // all vertex shaders define a main() function
         {
            gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
               // this line transforms the predefined attribute
               // gl_Vertex of type vec4 with the predefined
               // uniform gl_ModelViewProjectionMatrix of type mat4
               // and stores the result in the predefined output
               // variable gl_Position of type vec4.

            texCoord = gl_MultiTexCoord0;

         }

         #endif // here ends the definition of the vertex shader


         #ifdef FRAGMENT // here begins the fragment shader


         in vec4 texCoord;
         void main() // all fragment shaders define a main() function
         {
            vec3 in_rgb = texture(_ToneImg, texCoord.st).rgb*_MaxLum;

            vec3 out_rgb = clamp(mat3x3(_M_rec709_rgbdisp) * (in_rgb-_BlRGB), 0., 10000.);

            vec3 out_color;
            out_color.r = pow(out_rgb.r,_InvGamma.r);
            out_color.g = pow(out_rgb.g,_InvGamma.g);
            out_color.b = pow(out_rgb.b,_InvGamma.b);
            vec3 rgb_uc = out_color;

            int pattern_bit_x = int(gl_FragCoord.x) & 1;
            int pattern_bit_y = int(gl_FragCoord.y) & 1;
            int spatial_dither = (pattern_bit_x ^ pattern_bit_y);

            int bit_h = int(floor(out_color.g*_HWBitDepth*2.)) & 1;
            int bit_l = int(floor(out_color.g*_HWBitDepth*4.)) & 1;

            out_color = floor(out_color*_HWBitDepth)/_HWBitDepth;

            int full_bit = _FlickerCount & bit_h & bit_l;
            int half_bit = ((_FlickerCount+1) & bit_h) | (_FlickerCount & (bit_h^bit_l));
            out_color += vec3(1.,1.,1.)/_HWBitDepth * float(full_bit) +
                        vec3(1.,1.,1.)*float(spatial_dither)/_HWBitDepth * float(half_bit);

            if( _GamutWarning > 0 ) {

            if( (out_rgb.r < 0.) || (out_rgb.g < 0.) || (out_rgb.b < 0.) ) {
            out_color = vec3( 1., 0., 1. ) * float((_FlickerCount & 32)/32);
            }


            if( (out_rgb.r > 1.) || (out_rgb.g > 1.) || (out_rgb.b >= 1.)) {
            out_color = vec3( 1., 0., 0. ) * float((_FlickerCount & 32)/32);
            }

            }

    gl_FragColor.rgb = out_color;
    gl_FragColor.a = 1.;

         }


         #endif // here ends the definition of the fragment shader

         ENDGLSL // here ends the part in GLSL
      }

    }
}
