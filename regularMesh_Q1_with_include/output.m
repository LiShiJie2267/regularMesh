x=0:1/50:1;
y=0:1/40:1;
[X,Y]=meshgrid(x,y);
U=[0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , ;
0 , 0.00982922 , 0.0196197 , 0.0293327 , 0.0389299 , 0.0483735 , 0.0576262 , 0.0666515 , 0.0754137 , 0.0838783 , 0.0920119 , 0.0997823 , 0.107159 , 0.114113 , 0.120616 , 0.126643 , 0.132171 , 0.137177 , 0.141642 , 0.145547 , 0.148878 , 0.151622 , 0.153767 , 0.155306 , 0.156231 , 0.15654 , 0.156231 , 0.155306 , 0.153767 , 0.151622 , 0.148878 , 0.145547 , 0.141642 , 0.137177 , 0.132171 , 0.126643 , 0.120616 , 0.114113 , 0.107159 , 0.0997823 , 0.0920119 , 0.0838783 , 0.0754137 , 0.0666515 , 0.0576262 , 0.0483735 , 0.0389299 , 0.0293327 , 0.0196197 , 0.00982922 , 1.91577e-17 , ;
0 , 0.0194164 , 0.0387562 , 0.0579431 , 0.0769012 , 0.0955559 , 0.113833 , 0.131662 , 0.14897 , 0.165691 , 0.181758 , 0.197108 , 0.211679 , 0.225416 , 0.238262 , 0.250169 , 0.261088 , 0.270976 , 0.279795 , 0.28751 , 0.294091 , 0.29951 , 0.303748 , 0.306787 , 0.308615 , 0.309225 , 0.308615 , 0.306787 , 0.303748 , 0.29951 , 0.294091 , 0.28751 , 0.279795 , 0.270976 , 0.261088 , 0.250169 , 0.238262 , 0.225416 , 0.211679 , 0.197108 , 0.181758 , 0.165691 , 0.14897 , 0.131662 , 0.113833 , 0.0955559 , 0.0769012 , 0.0579431 , 0.0387562 , 0.0194164 , 3.78437e-17 , ;
0 , 0.0285255 , 0.0569385 , 0.0851267 , 0.112979 , 0.140385 , 0.167238 , 0.19343 , 0.218859 , 0.243424 , 0.267029 , 0.28958 , 0.310987 , 0.331168 , 0.350042 , 0.367534 , 0.383575 , 0.398103 , 0.41106 , 0.422394 , 0.432062 , 0.440024 , 0.44625 , 0.450714 , 0.4534 , 0.454297 , 0.4534 , 0.450714 , 0.44625 , 0.440024 , 0.432062 , 0.422394 , 0.41106 , 0.398103 , 0.383575 , 0.367534 , 0.350042 , 0.331168 , 0.310987 , 0.28958 , 0.267029 , 0.243424 , 0.218859 , 0.19343 , 0.167238 , 0.140385 , 0.112979 , 0.0851267 , 0.0569385 , 0.0285255 , 5.55978e-17 , ;
0 , 0.0369322 , 0.0737187 , 0.110214 , 0.146275 , 0.181758 , 0.216524 , 0.250436 , 0.283359 , 0.315163 , 0.345724 , 0.374921 , 0.402638 , 0.428766 , 0.453202 , 0.475849 , 0.496618 , 0.515427 , 0.532203 , 0.546877 , 0.559394 , 0.569703 , 0.577763 , 0.583544 , 0.587021 , 0.588182 , 0.587021 , 0.583544 , 0.577763 , 0.569703 , 0.559394 , 0.546877 , 0.532203 , 0.515427 , 0.496618 , 0.475849 , 0.453202 , 0.428766 , 0.402638 , 0.374921 , 0.345724 , 0.315163 , 0.283359 , 0.250436 , 0.216524 , 0.181758 , 0.146275 , 0.110214 , 0.0737187 , 0.0369322 , 7.19829e-17 , ;
0 , 0.0444295 , 0.0886837 , 0.132588 , 0.175969 , 0.218655 , 0.260479 , 0.301274 , 0.340881 , 0.379142 , 0.415907 , 0.451031 , 0.484374 , 0.515806 , 0.545203 , 0.572447 , 0.597433 , 0.62006 , 0.640241 , 0.657895 , 0.672952 , 0.685354 , 0.69505 , 0.702004 , 0.706187 , 0.707584 , 0.706187 , 0.702004 , 0.69505 , 0.685354 , 0.672952 , 0.657895 , 0.640241 , 0.62006 , 0.597433 , 0.572447 , 0.545203 , 0.515806 , 0.484374 , 0.451031 , 0.415907 , 0.379142 , 0.340881 , 0.301274 , 0.260479 , 0.218655 , 0.175969 , 0.132588 , 0.0886837 , 0.0444295 , 8.65956e-17 , ;
0 , 0.0508329 , 0.101465 , 0.151697 , 0.20133 , 0.250169 , 0.29802 , 0.344695 , 0.39001 , 0.433785 , 0.475849 , 0.516035 , 0.554184 , 0.590146 , 0.623779 , 0.65495 , 0.683536 , 0.709425 , 0.732514 , 0.752712 , 0.76994 , 0.784129 , 0.795223 , 0.803179 , 0.807965 , 0.809563 , 0.807965 , 0.803179 , 0.795223 , 0.784129 , 0.76994 , 0.752712 , 0.732514 , 0.709425 , 0.683536 , 0.65495 , 0.623779 , 0.590146 , 0.554184 , 0.516035 , 0.475849 , 0.433785 , 0.39001 , 0.344695 , 0.29802 , 0.250169 , 0.20133 , 0.151697 , 0.101465 , 0.0508329 , 9.9076e-17 , ;
0 , 0.0559845 , 0.111748 , 0.167071 , 0.221734 , 0.275522 , 0.328223 , 0.379628 , 0.429535 , 0.477747 , 0.524074 , 0.568332 , 0.610347 , 0.649954 , 0.686995 , 0.721325 , 0.752809 , 0.781321 , 0.80675 , 0.828996 , 0.847969 , 0.863596 , 0.875815 , 0.884577 , 0.889848 , 0.891607 , 0.889848 , 0.884577 , 0.875815 , 0.863596 , 0.847969 , 0.828996 , 0.80675 , 0.781321 , 0.752809 , 0.721325 , 0.686995 , 0.649954 , 0.610347 , 0.568332 , 0.524074 , 0.477747 , 0.429535 , 0.379628 , 0.328223 , 0.275522 , 0.221734 , 0.167071 , 0.111748 , 0.0559845 , 1.09117e-16 , ;
0 , 0.0597576 , 0.119279 , 0.17833 , 0.236678 , 0.294091 , 0.350343 , 0.405213 , 0.458484 , 0.509945 , 0.559394 , 0.606635 , 0.651482 , 0.693758 , 0.733296 , 0.76994 , 0.803545 , 0.833979 , 0.861122 , 0.884866 , 0.905118 , 0.921798 , 0.934841 , 0.944193 , 0.94982 , 0.951698 , 0.94982 , 0.944193 , 0.934841 , 0.921798 , 0.905118 , 0.884866 , 0.861122 , 0.833979 , 0.803545 , 0.76994 , 0.733296 , 0.693758 , 0.651482 , 0.606635 , 0.559394 , 0.509945 , 0.458484 , 0.405213 , 0.350343 , 0.294091 , 0.236678 , 0.17833 , 0.119279 , 0.0597576 , 1.16471e-16 , ;
0 , 0.0620593 , 0.123874 , 0.185199 , 0.245794 , 0.305418 , 0.363838 , 0.420821 , 0.476143 , 0.529587 , 0.58094 , 0.630001 , 0.676575 , 0.720479 , 0.76154 , 0.799595 , 0.834495 , 0.866102 , 0.89429 , 0.918949 , 0.939981 , 0.957303 , 0.970848 , 0.980561 , 0.986404 , 0.988354 , 0.986404 , 0.980561 , 0.970848 , 0.957303 , 0.939981 , 0.918949 , 0.89429 , 0.866102 , 0.834495 , 0.799595 , 0.76154 , 0.720479 , 0.676575 , 0.630001 , 0.58094 , 0.529587 , 0.476143 , 0.420821 , 0.363838 , 0.305418 , 0.245794 , 0.185199 , 0.123874 , 0.0620593 , 1.20957e-16 , ;
0 , 0.0628329 , 0.125418 , 0.187508 , 0.248858 , 0.309225 , 0.368373 , 0.426066 , 0.482079 , 0.536188 , 0.588182 , 0.637854 , 0.685009 , 0.72946 , 0.771033 , 0.809563 , 0.844897 , 0.876898 , 0.905437 , 0.930403 , 0.951698 , 0.969236 , 0.98295 , 0.992784 , 0.9987 , 1.00067 , 0.9987 , 0.992784 , 0.98295 , 0.969236 , 0.951698 , 0.930403 , 0.905437 , 0.876898 , 0.844897 , 0.809563 , 0.771033 , 0.72946 , 0.685009 , 0.637854 , 0.588182 , 0.536188 , 0.482079 , 0.426066 , 0.368373 , 0.309225 , 0.248858 , 0.187508 , 0.125418 , 0.0628329 , 1.22465e-16 , ;
0 , 0.0620593 , 0.123874 , 0.185199 , 0.245794 , 0.305418 , 0.363838 , 0.420821 , 0.476143 , 0.529587 , 0.58094 , 0.630001 , 0.676575 , 0.720479 , 0.76154 , 0.799595 , 0.834495 , 0.866102 , 0.89429 , 0.918949 , 0.939981 , 0.957303 , 0.970848 , 0.980561 , 0.986404 , 0.988354 , 0.986404 , 0.980561 , 0.970848 , 0.957303 , 0.939981 , 0.918949 , 0.89429 , 0.866102 , 0.834495 , 0.799595 , 0.76154 , 0.720479 , 0.676575 , 0.630001 , 0.58094 , 0.529587 , 0.476143 , 0.420821 , 0.363838 , 0.305418 , 0.245794 , 0.185199 , 0.123874 , 0.0620593 , 1.20957e-16 , ;
0 , 0.0597576 , 0.119279 , 0.17833 , 0.236678 , 0.294091 , 0.350343 , 0.405213 , 0.458484 , 0.509945 , 0.559394 , 0.606635 , 0.651482 , 0.693758 , 0.733296 , 0.76994 , 0.803545 , 0.833979 , 0.861122 , 0.884866 , 0.905118 , 0.921798 , 0.934841 , 0.944193 , 0.94982 , 0.951698 , 0.94982 , 0.944193 , 0.934841 , 0.921798 , 0.905118 , 0.884866 , 0.861122 , 0.833979 , 0.803545 , 0.76994 , 0.733296 , 0.693758 , 0.651482 , 0.606635 , 0.559394 , 0.509945 , 0.458484 , 0.405213 , 0.350343 , 0.294091 , 0.236678 , 0.17833 , 0.119279 , 0.0597576 , 1.16471e-16 , ;
0 , 0.0559845 , 0.111748 , 0.167071 , 0.221734 , 0.275522 , 0.328223 , 0.379628 , 0.429535 , 0.477747 , 0.524074 , 0.568332 , 0.610347 , 0.649954 , 0.686995 , 0.721325 , 0.752809 , 0.781321 , 0.80675 , 0.828996 , 0.847969 , 0.863596 , 0.875815 , 0.884577 , 0.889848 , 0.891607 , 0.889848 , 0.884577 , 0.875815 , 0.863596 , 0.847969 , 0.828996 , 0.80675 , 0.781321 , 0.752809 , 0.721325 , 0.686995 , 0.649954 , 0.610347 , 0.568332 , 0.524074 , 0.477747 , 0.429535 , 0.379628 , 0.328223 , 0.275522 , 0.221734 , 0.167071 , 0.111748 , 0.0559845 , 1.09117e-16 , ;
0 , 0.0508329 , 0.101465 , 0.151697 , 0.20133 , 0.250169 , 0.29802 , 0.344695 , 0.39001 , 0.433785 , 0.475849 , 0.516035 , 0.554184 , 0.590146 , 0.623779 , 0.65495 , 0.683536 , 0.709425 , 0.732514 , 0.752712 , 0.76994 , 0.784129 , 0.795223 , 0.803179 , 0.807965 , 0.809563 , 0.807965 , 0.803179 , 0.795223 , 0.784129 , 0.76994 , 0.752712 , 0.732514 , 0.709425 , 0.683536 , 0.65495 , 0.623779 , 0.590146 , 0.554184 , 0.516035 , 0.475849 , 0.433785 , 0.39001 , 0.344695 , 0.29802 , 0.250169 , 0.20133 , 0.151697 , 0.101465 , 0.0508329 , 9.9076e-17 , ;
0 , 0.0444295 , 0.0886837 , 0.132588 , 0.175969 , 0.218655 , 0.260479 , 0.301274 , 0.340881 , 0.379142 , 0.415907 , 0.451031 , 0.484374 , 0.515806 , 0.545203 , 0.572447 , 0.597433 , 0.62006 , 0.640241 , 0.657895 , 0.672952 , 0.685354 , 0.69505 , 0.702004 , 0.706187 , 0.707584 , 0.706187 , 0.702004 , 0.69505 , 0.685354 , 0.672952 , 0.657895 , 0.640241 , 0.62006 , 0.597433 , 0.572447 , 0.545203 , 0.515806 , 0.484374 , 0.451031 , 0.415907 , 0.379142 , 0.340881 , 0.301274 , 0.260479 , 0.218655 , 0.175969 , 0.132588 , 0.0886837 , 0.0444295 , 8.65956e-17 , ;
0 , 0.0369322 , 0.0737187 , 0.110214 , 0.146275 , 0.181758 , 0.216524 , 0.250436 , 0.283359 , 0.315163 , 0.345724 , 0.374921 , 0.402638 , 0.428766 , 0.453202 , 0.475849 , 0.496618 , 0.515427 , 0.532203 , 0.546877 , 0.559394 , 0.569703 , 0.577763 , 0.583544 , 0.587021 , 0.588182 , 0.587021 , 0.583544 , 0.577763 , 0.569703 , 0.559394 , 0.546877 , 0.532203 , 0.515427 , 0.496618 , 0.475849 , 0.453202 , 0.428766 , 0.402638 , 0.374921 , 0.345724 , 0.315163 , 0.283359 , 0.250436 , 0.216524 , 0.181758 , 0.146275 , 0.110214 , 0.0737187 , 0.0369322 , 7.19829e-17 , ;
0 , 0.0285255 , 0.0569385 , 0.0851267 , 0.112979 , 0.140385 , 0.167238 , 0.19343 , 0.218859 , 0.243424 , 0.267029 , 0.28958 , 0.310987 , 0.331168 , 0.350042 , 0.367534 , 0.383575 , 0.398103 , 0.41106 , 0.422394 , 0.432062 , 0.440024 , 0.44625 , 0.450714 , 0.4534 , 0.454297 , 0.4534 , 0.450714 , 0.44625 , 0.440024 , 0.432062 , 0.422394 , 0.41106 , 0.398103 , 0.383575 , 0.367534 , 0.350042 , 0.331168 , 0.310987 , 0.28958 , 0.267029 , 0.243424 , 0.218859 , 0.19343 , 0.167238 , 0.140385 , 0.112979 , 0.0851267 , 0.0569385 , 0.0285255 , 5.55978e-17 , ;
0 , 0.0194164 , 0.0387562 , 0.0579431 , 0.0769012 , 0.0955559 , 0.113833 , 0.131662 , 0.14897 , 0.165691 , 0.181758 , 0.197108 , 0.211679 , 0.225416 , 0.238262 , 0.250169 , 0.261088 , 0.270976 , 0.279795 , 0.28751 , 0.294091 , 0.29951 , 0.303748 , 0.306787 , 0.308615 , 0.309225 , 0.308615 , 0.306787 , 0.303748 , 0.29951 , 0.294091 , 0.28751 , 0.279795 , 0.270976 , 0.261088 , 0.250169 , 0.238262 , 0.225416 , 0.211679 , 0.197108 , 0.181758 , 0.165691 , 0.14897 , 0.131662 , 0.113833 , 0.0955559 , 0.0769012 , 0.0579431 , 0.0387562 , 0.0194164 , 3.78437e-17 , ;
0 , 0.00982922 , 0.0196197 , 0.0293327 , 0.0389299 , 0.0483735 , 0.0576262 , 0.0666515 , 0.0754137 , 0.0838783 , 0.0920119 , 0.0997823 , 0.107159 , 0.114113 , 0.120616 , 0.126643 , 0.132171 , 0.137177 , 0.141642 , 0.145547 , 0.148878 , 0.151622 , 0.153767 , 0.155306 , 0.156231 , 0.15654 , 0.156231 , 0.155306 , 0.153767 , 0.151622 , 0.148878 , 0.145547 , 0.141642 , 0.137177 , 0.132171 , 0.126643 , 0.120616 , 0.114113 , 0.107159 , 0.0997823 , 0.0920119 , 0.0838783 , 0.0754137 , 0.0666515 , 0.0576262 , 0.0483735 , 0.0389299 , 0.0293327 , 0.0196197 , 0.00982922 , 1.91577e-17 , ;
0 , -2.23148e-14 , -4.44348e-14 , -6.57021e-14 , -8.59374e-14 , -1.04266e-13 , -1.20607e-13 , -1.34279e-13 , -1.45585e-13 , -1.54224e-13 , -1.60794e-13 , -1.65241e-13 , -1.68295e-13 , -1.69907e-13 , -1.70705e-13 , -1.70585e-13 , -1.70052e-13 , -1.68988e-13 , -1.67724e-13 , -1.6611e-13 , -1.64416e-13 , -1.62477e-13 , -1.60474e-13 , -1.58265e-13 , -1.55998e-13 , -1.53536e-13 , -1.50992e-13 , -1.48228e-13 , -1.45387e-13 , -1.42254e-13 , -1.39048e-13 , -1.35469e-13 , -1.31752e-13 , -1.27515e-13 , -1.23034e-13 , -1.17856e-13 , -1.12245e-13 , -1.05807e-13 , -9.88108e-14 , -9.10553e-14 , -8.2768e-14 , -7.40212e-14 , -6.49659e-14 , -5.58979e-14 , -4.68896e-14 , -3.82532e-14 , -2.99418e-14 , -2.20816e-14 , -1.45357e-14 , -7.22778e-15 , 1.49976e-32 , ;
0 , -0.00982922 , -0.0196197 , -0.0293327 , -0.0389299 , -0.0483735 , -0.0576262 , -0.0666515 , -0.0754137 , -0.0838783 , -0.0920119 , -0.0997823 , -0.107159 , -0.114113 , -0.120616 , -0.126643 , -0.132171 , -0.137177 , -0.141642 , -0.145547 , -0.148878 , -0.151622 , -0.153767 , -0.155306 , -0.156231 , -0.15654 , -0.156231 , -0.155306 , -0.153767 , -0.151622 , -0.148878 , -0.145547 , -0.141642 , -0.137177 , -0.132171 , -0.126643 , -0.120616 , -0.114113 , -0.107159 , -0.0997823 , -0.0920119 , -0.0838783 , -0.0754137 , -0.0666515 , -0.0576262 , -0.0483735 , -0.0389299 , -0.0293327 , -0.0196197 , -0.00982922 , -1.91577e-17 , ;
0 , -0.0194164 , -0.0387562 , -0.0579431 , -0.0769012 , -0.0955559 , -0.113833 , -0.131662 , -0.14897 , -0.165691 , -0.181758 , -0.197108 , -0.211679 , -0.225416 , -0.238262 , -0.250169 , -0.261088 , -0.270976 , -0.279795 , -0.28751 , -0.294091 , -0.29951 , -0.303748 , -0.306787 , -0.308615 , -0.309225 , -0.308615 , -0.306787 , -0.303748 , -0.29951 , -0.294091 , -0.28751 , -0.279795 , -0.270976 , -0.261088 , -0.250169 , -0.238262 , -0.225416 , -0.211679 , -0.197108 , -0.181758 , -0.165691 , -0.14897 , -0.131662 , -0.113833 , -0.0955559 , -0.0769012 , -0.0579431 , -0.0387562 , -0.0194164 , -3.78437e-17 , ;
0 , -0.0285255 , -0.0569385 , -0.0851267 , -0.112979 , -0.140385 , -0.167238 , -0.19343 , -0.218859 , -0.243424 , -0.267029 , -0.28958 , -0.310987 , -0.331168 , -0.350042 , -0.367534 , -0.383575 , -0.398103 , -0.41106 , -0.422394 , -0.432062 , -0.440024 , -0.44625 , -0.450714 , -0.4534 , -0.454297 , -0.4534 , -0.450714 , -0.44625 , -0.440024 , -0.432062 , -0.422394 , -0.41106 , -0.398103 , -0.383575 , -0.367534 , -0.350042 , -0.331168 , -0.310987 , -0.28958 , -0.267029 , -0.243424 , -0.218859 , -0.19343 , -0.167238 , -0.140385 , -0.112979 , -0.0851267 , -0.0569385 , -0.0285255 , -5.55978e-17 , ;
0 , -0.0369322 , -0.0737187 , -0.110214 , -0.146275 , -0.181758 , -0.216524 , -0.250436 , -0.283359 , -0.315163 , -0.345724 , -0.374921 , -0.402638 , -0.428766 , -0.453202 , -0.475849 , -0.496618 , -0.515427 , -0.532203 , -0.546877 , -0.559394 , -0.569703 , -0.577763 , -0.583544 , -0.587021 , -0.588182 , -0.587021 , -0.583544 , -0.577763 , -0.569703 , -0.559394 , -0.546877 , -0.532203 , -0.515427 , -0.496618 , -0.475849 , -0.453202 , -0.428766 , -0.402638 , -0.374921 , -0.345724 , -0.315163 , -0.283359 , -0.250436 , -0.216524 , -0.181758 , -0.146275 , -0.110214 , -0.0737187 , -0.0369322 , -7.19829e-17 , ;
0 , -0.0444295 , -0.0886837 , -0.132588 , -0.175969 , -0.218655 , -0.260479 , -0.301274 , -0.340881 , -0.379142 , -0.415907 , -0.451031 , -0.484374 , -0.515806 , -0.545203 , -0.572447 , -0.597433 , -0.62006 , -0.640241 , -0.657895 , -0.672952 , -0.685354 , -0.69505 , -0.702004 , -0.706187 , -0.707584 , -0.706187 , -0.702004 , -0.69505 , -0.685354 , -0.672952 , -0.657895 , -0.640241 , -0.62006 , -0.597433 , -0.572447 , -0.545203 , -0.515806 , -0.484374 , -0.451031 , -0.415907 , -0.379142 , -0.340881 , -0.301274 , -0.260479 , -0.218655 , -0.175969 , -0.132588 , -0.0886837 , -0.0444295 , -8.65956e-17 , ;
0 , -0.0508329 , -0.101465 , -0.151697 , -0.20133 , -0.250169 , -0.29802 , -0.344695 , -0.39001 , -0.433785 , -0.475849 , -0.516035 , -0.554184 , -0.590146 , -0.623779 , -0.65495 , -0.683536 , -0.709425 , -0.732514 , -0.752712 , -0.76994 , -0.784129 , -0.795223 , -0.803179 , -0.807965 , -0.809563 , -0.807965 , -0.803179 , -0.795223 , -0.784129 , -0.76994 , -0.752712 , -0.732514 , -0.709425 , -0.683536 , -0.65495 , -0.623779 , -0.590146 , -0.554184 , -0.516035 , -0.475849 , -0.433785 , -0.39001 , -0.344695 , -0.29802 , -0.250169 , -0.20133 , -0.151697 , -0.101465 , -0.0508329 , -9.9076e-17 , ;
0 , -0.0559845 , -0.111748 , -0.167071 , -0.221734 , -0.275522 , -0.328223 , -0.379628 , -0.429535 , -0.477747 , -0.524074 , -0.568332 , -0.610347 , -0.649954 , -0.686995 , -0.721325 , -0.752809 , -0.781321 , -0.80675 , -0.828996 , -0.847969 , -0.863596 , -0.875815 , -0.884577 , -0.889848 , -0.891607 , -0.889848 , -0.884577 , -0.875815 , -0.863596 , -0.847969 , -0.828996 , -0.80675 , -0.781321 , -0.752809 , -0.721325 , -0.686995 , -0.649954 , -0.610347 , -0.568332 , -0.524074 , -0.477747 , -0.429535 , -0.379628 , -0.328223 , -0.275522 , -0.221734 , -0.167071 , -0.111748 , -0.0559845 , -1.09117e-16 , ;
0 , -0.0597576 , -0.119279 , -0.17833 , -0.236678 , -0.294091 , -0.350343 , -0.405213 , -0.458484 , -0.509945 , -0.559394 , -0.606635 , -0.651482 , -0.693758 , -0.733296 , -0.76994 , -0.803545 , -0.833979 , -0.861122 , -0.884866 , -0.905118 , -0.921798 , -0.934841 , -0.944193 , -0.94982 , -0.951698 , -0.94982 , -0.944193 , -0.934841 , -0.921798 , -0.905118 , -0.884866 , -0.861122 , -0.833979 , -0.803545 , -0.76994 , -0.733296 , -0.693758 , -0.651482 , -0.606635 , -0.559394 , -0.509945 , -0.458484 , -0.405213 , -0.350343 , -0.294091 , -0.236678 , -0.17833 , -0.119279 , -0.0597576 , -1.16471e-16 , ;
0 , -0.0620593 , -0.123874 , -0.185199 , -0.245794 , -0.305418 , -0.363838 , -0.420821 , -0.476143 , -0.529587 , -0.58094 , -0.630001 , -0.676575 , -0.720479 , -0.76154 , -0.799595 , -0.834495 , -0.866102 , -0.89429 , -0.918949 , -0.939981 , -0.957303 , -0.970848 , -0.980561 , -0.986404 , -0.988354 , -0.986404 , -0.980561 , -0.970848 , -0.957303 , -0.939981 , -0.918949 , -0.89429 , -0.866102 , -0.834495 , -0.799595 , -0.76154 , -0.720479 , -0.676575 , -0.630001 , -0.58094 , -0.529587 , -0.476143 , -0.420821 , -0.363838 , -0.305418 , -0.245794 , -0.185199 , -0.123874 , -0.0620593 , -1.20957e-16 , ;
0 , -0.0628329 , -0.125418 , -0.187508 , -0.248858 , -0.309225 , -0.368373 , -0.426066 , -0.482079 , -0.536188 , -0.588182 , -0.637854 , -0.685009 , -0.72946 , -0.771033 , -0.809563 , -0.844897 , -0.876898 , -0.905437 , -0.930403 , -0.951698 , -0.969236 , -0.98295 , -0.992784 , -0.9987 , -1.00067 , -0.9987 , -0.992784 , -0.98295 , -0.969236 , -0.951698 , -0.930403 , -0.905437 , -0.876898 , -0.844897 , -0.809563 , -0.771033 , -0.72946 , -0.685009 , -0.637854 , -0.588182 , -0.536188 , -0.482079 , -0.426066 , -0.368373 , -0.309225 , -0.248858 , -0.187508 , -0.125418 , -0.0628329 , -1.22465e-16 , ;
0 , -0.0620593 , -0.123874 , -0.185199 , -0.245794 , -0.305418 , -0.363838 , -0.420821 , -0.476143 , -0.529587 , -0.58094 , -0.630001 , -0.676575 , -0.720479 , -0.76154 , -0.799595 , -0.834495 , -0.866102 , -0.89429 , -0.918949 , -0.939981 , -0.957303 , -0.970848 , -0.980561 , -0.986404 , -0.988354 , -0.986404 , -0.980561 , -0.970848 , -0.957303 , -0.939981 , -0.918949 , -0.89429 , -0.866102 , -0.834495 , -0.799595 , -0.76154 , -0.720479 , -0.676575 , -0.630001 , -0.58094 , -0.529587 , -0.476143 , -0.420821 , -0.363838 , -0.305418 , -0.245794 , -0.185199 , -0.123874 , -0.0620593 , -1.20957e-16 , ;
0 , -0.0597576 , -0.119279 , -0.17833 , -0.236678 , -0.294091 , -0.350343 , -0.405213 , -0.458484 , -0.509945 , -0.559394 , -0.606635 , -0.651482 , -0.693758 , -0.733296 , -0.76994 , -0.803545 , -0.833979 , -0.861122 , -0.884866 , -0.905118 , -0.921798 , -0.934841 , -0.944193 , -0.94982 , -0.951698 , -0.94982 , -0.944193 , -0.934841 , -0.921798 , -0.905118 , -0.884866 , -0.861122 , -0.833979 , -0.803545 , -0.76994 , -0.733296 , -0.693758 , -0.651482 , -0.606635 , -0.559394 , -0.509945 , -0.458484 , -0.405213 , -0.350343 , -0.294091 , -0.236678 , -0.17833 , -0.119279 , -0.0597576 , -1.16471e-16 , ;
0 , -0.0559845 , -0.111748 , -0.167071 , -0.221734 , -0.275522 , -0.328223 , -0.379628 , -0.429535 , -0.477747 , -0.524074 , -0.568332 , -0.610347 , -0.649954 , -0.686995 , -0.721325 , -0.752809 , -0.781321 , -0.80675 , -0.828996 , -0.847969 , -0.863596 , -0.875815 , -0.884577 , -0.889848 , -0.891607 , -0.889848 , -0.884577 , -0.875815 , -0.863596 , -0.847969 , -0.828996 , -0.80675 , -0.781321 , -0.752809 , -0.721325 , -0.686995 , -0.649954 , -0.610347 , -0.568332 , -0.524074 , -0.477747 , -0.429535 , -0.379628 , -0.328223 , -0.275522 , -0.221734 , -0.167071 , -0.111748 , -0.0559845 , -1.09117e-16 , ;
0 , -0.0508329 , -0.101465 , -0.151697 , -0.20133 , -0.250169 , -0.29802 , -0.344695 , -0.39001 , -0.433785 , -0.475849 , -0.516035 , -0.554184 , -0.590146 , -0.623779 , -0.65495 , -0.683536 , -0.709425 , -0.732514 , -0.752712 , -0.76994 , -0.784129 , -0.795223 , -0.803179 , -0.807965 , -0.809563 , -0.807965 , -0.803179 , -0.795223 , -0.784129 , -0.76994 , -0.752712 , -0.732514 , -0.709425 , -0.683536 , -0.65495 , -0.623779 , -0.590146 , -0.554184 , -0.516035 , -0.475849 , -0.433785 , -0.39001 , -0.344695 , -0.29802 , -0.250169 , -0.20133 , -0.151697 , -0.101465 , -0.0508329 , -9.9076e-17 , ;
0 , -0.0444295 , -0.0886837 , -0.132588 , -0.175969 , -0.218655 , -0.260479 , -0.301274 , -0.340881 , -0.379142 , -0.415907 , -0.451031 , -0.484374 , -0.515806 , -0.545203 , -0.572447 , -0.597433 , -0.62006 , -0.640241 , -0.657895 , -0.672952 , -0.685354 , -0.69505 , -0.702004 , -0.706187 , -0.707584 , -0.706187 , -0.702004 , -0.69505 , -0.685354 , -0.672952 , -0.657895 , -0.640241 , -0.62006 , -0.597433 , -0.572447 , -0.545203 , -0.515806 , -0.484374 , -0.451031 , -0.415907 , -0.379142 , -0.340881 , -0.301274 , -0.260479 , -0.218655 , -0.175969 , -0.132588 , -0.0886837 , -0.0444295 , -8.65956e-17 , ;
0 , -0.0369322 , -0.0737187 , -0.110214 , -0.146275 , -0.181758 , -0.216524 , -0.250436 , -0.283359 , -0.315163 , -0.345724 , -0.374921 , -0.402638 , -0.428766 , -0.453202 , -0.475849 , -0.496618 , -0.515427 , -0.532203 , -0.546877 , -0.559394 , -0.569703 , -0.577763 , -0.583544 , -0.587021 , -0.588182 , -0.587021 , -0.583544 , -0.577763 , -0.569703 , -0.559394 , -0.546877 , -0.532203 , -0.515427 , -0.496618 , -0.475849 , -0.453202 , -0.428766 , -0.402638 , -0.374921 , -0.345724 , -0.315163 , -0.283359 , -0.250436 , -0.216524 , -0.181758 , -0.146275 , -0.110214 , -0.0737187 , -0.0369322 , -7.19829e-17 , ;
0 , -0.0285255 , -0.0569385 , -0.0851267 , -0.112979 , -0.140385 , -0.167238 , -0.19343 , -0.218859 , -0.243424 , -0.267029 , -0.28958 , -0.310987 , -0.331168 , -0.350042 , -0.367534 , -0.383575 , -0.398103 , -0.41106 , -0.422394 , -0.432062 , -0.440024 , -0.44625 , -0.450714 , -0.4534 , -0.454297 , -0.4534 , -0.450714 , -0.44625 , -0.440024 , -0.432062 , -0.422394 , -0.41106 , -0.398103 , -0.383575 , -0.367534 , -0.350042 , -0.331168 , -0.310987 , -0.28958 , -0.267029 , -0.243424 , -0.218859 , -0.19343 , -0.167238 , -0.140385 , -0.112979 , -0.0851267 , -0.0569385 , -0.0285255 , -5.55978e-17 , ;
0 , -0.0194164 , -0.0387562 , -0.0579431 , -0.0769012 , -0.0955559 , -0.113833 , -0.131662 , -0.14897 , -0.165691 , -0.181758 , -0.197108 , -0.211679 , -0.225416 , -0.238262 , -0.250169 , -0.261088 , -0.270976 , -0.279795 , -0.28751 , -0.294091 , -0.29951 , -0.303748 , -0.306787 , -0.308615 , -0.309225 , -0.308615 , -0.306787 , -0.303748 , -0.29951 , -0.294091 , -0.28751 , -0.279795 , -0.270976 , -0.261088 , -0.250169 , -0.238262 , -0.225416 , -0.211679 , -0.197108 , -0.181758 , -0.165691 , -0.14897 , -0.131662 , -0.113833 , -0.0955559 , -0.0769012 , -0.0579431 , -0.0387562 , -0.0194164 , -3.78437e-17 , ;
0 , -0.00982922 , -0.0196197 , -0.0293327 , -0.0389299 , -0.0483735 , -0.0576262 , -0.0666515 , -0.0754137 , -0.0838783 , -0.0920119 , -0.0997823 , -0.107159 , -0.114113 , -0.120616 , -0.126643 , -0.132171 , -0.137177 , -0.141642 , -0.145547 , -0.148878 , -0.151622 , -0.153767 , -0.155306 , -0.156231 , -0.15654 , -0.156231 , -0.155306 , -0.153767 , -0.151622 , -0.148878 , -0.145547 , -0.141642 , -0.137177 , -0.132171 , -0.126643 , -0.120616 , -0.114113 , -0.107159 , -0.0997823 , -0.0920119 , -0.0838783 , -0.0754137 , -0.0666515 , -0.0576262 , -0.0483735 , -0.0389299 , -0.0293327 , -0.0196197 , -0.00982922 , -1.91577e-17 , ;
0 , -1.53792e-17 , -3.06978e-17 , -4.58952e-17 , -6.09115e-17 , -7.56873e-17 , -9.01645e-17 , -1.04286e-16 , -1.17996e-16 , -1.3124e-16 , -1.43966e-16 , -1.56124e-16 , -1.67666e-16 , -1.78546e-16 , -1.88721e-16 , -1.98152e-16 , -2.06801e-16 , -2.14633e-16 , -2.21619e-16 , -2.2773e-16 , -2.32942e-16 , -2.37234e-16 , -2.40591e-16 , -2.42998e-16 , -2.44446e-16 , -2.44929e-16 , -2.44446e-16 , -2.42998e-16 , -2.40591e-16 , -2.37234e-16 , -2.32942e-16 , -2.2773e-16 , -2.21619e-16 , -2.14633e-16 , -2.06801e-16 , -1.98152e-16 , -1.88721e-16 , -1.78546e-16 , -1.67666e-16 , -1.56124e-16 , -1.43966e-16 , -1.3124e-16 , -1.17996e-16 , -1.04286e-16 , -9.01645e-17 , -7.56873e-17 , -6.09115e-17 , -4.58952e-17 , -3.06978e-17 , -1.53792e-17 , -2.99952e-32 , ;
]
surf(x,y,U);
