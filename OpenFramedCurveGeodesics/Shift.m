function shift_vector = Shift(piece)
    eps=0.01;
    [~,n]=size(piece);
    shift_vector = [piece(1,n)*ones(1,100);piece(2,n)*ones(1,100);piece(3,n)*ones(1,100)]+eps*[ones(1,100);zeros(1,100);zeros(1,100)];