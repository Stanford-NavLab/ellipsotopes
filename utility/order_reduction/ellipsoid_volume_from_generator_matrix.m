function v = ellipsoid_volume_from_generator_matrix(G)
    Gi = pinv(G) ;
    v = det(inv(Gi'*Gi)) ;
end