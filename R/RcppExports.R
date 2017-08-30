
mh_sample_ne1_2  <- function ( fwdne, fwdnco , tau, prop_sigma, dh , lterms ){
	.Call( 'phyloma_mh_sample_ne1_2', PACKAGE = 'armasky', fwdne, fwdnco, tau, prop_sigma, dh, lterms )
}

mh_sample_ne2_2  <- function ( fwdne, fwdnco , tau, prop_sigma, dh, zxb, lterms ){
	.Call( 'phyloma_mh_sample_ne2_2', PACKAGE = 'armasky', fwdne, fwdnco, tau, prop_sigma, dh , zxb, lterms )
}
