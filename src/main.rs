// Program to simulate quaNTum wavepacket
use std::fs;
use std::io::Write;

// hbar^2/m0, where m0 is electron mass
const H2M: f64 = 0.0762; //eV * nm^2
const PIE: f64 = std::f64::consts::PI;
// effective mass
const ME: f64 = 0.067;

// wavepacket parameters
const LAM: f64 = 3.0; //nm
const X0: f64 = 20.0; //nm
const A: f64 = 3.0; //nm

// auxillary parameters
const K: f64 = 2.0*PIE/LAM; 
const A_A: f64 = H2M/ME/2.0;
const V: f64 = 2.0*A_A*K; // eV * nm

// grid parameters
const L_L: f64 = 120.0; //nm
const DX: f64 = 0.05; //nm
const NX: usize = (L_L/DX) as usize;
const DT: f64 = (DX*DX)/6.0/A_A; //1/eV
const T_T: f64 = L_L/V;
const NT: usize = (T_T/DT) as usize;
const ST: i32 = (NT as f64/200.0) as i32;
const B_B: f64 = 2.0*A_A*DT/(DX*DX);

// parameters of the potential
const UU0: f64 = 2.0; //eV
const U0: f64 = UU0*2.0*DT; // dimensionless
const XP: f64 = 45.0; //nm
const AP: f64 = 2.5; //nm
const TOL: f64 = 1.0e-8;

fn main() {
    //make directory
    let dirname = "results";
    fs::create_dir(&dirname).expect("Error creating directory");

    //initiate the x array
    let x: [f64; NX] = core::array::from_fn(|i| i as f64 * DX);

    fn ux(x: f64) -> f64 {
        let p = (-(x-XP)*(x-XP)/2.0/AP/AP).exp();
        if p > TOL {
            U0 * p
        } else { 
            0.0
        }
    }

    // Potential array
    let uj: [f64; NX] = core::array::from_fn(|i| ux(x[i]));


    // Initial state arrays
    let mut psi0: [f64; NX] = [0.0; NX];
    let mut phi0: [f64; NX] = [0.0; NX];
    let c_c = (0.5/PIE/A).sqrt();
    for j in 0..NX {
        let p = (-(x[j]-X0)*(x[j]-X0)/2.0/A/A).exp();
        if p > TOL{
            psi0[j] = c_c * p * (K * x[j]).cos();
            phi0[j] = c_c * p * (K * x[j]).sin();
        }
    }

    let mut psi1 = psi0.clone();
    let mut phi1 = phi0.clone();
    let mut psi2 = psi0.clone();
    let mut phi2 = phi0.clone();  

    // Euler step
    psi1[0] =  psi0[0] - 0.5*B_B*(- 2.0*phi0[0] + phi0[1]); //+ 0.5*uj[0]*phi0[0];
    phi1[0] =  phi0[0] + 0.5*B_B*(- 2.0*psi0[0] + psi0[1]); //- 0.5*uj[0]*psi0[0];   
    for j in 1..(NX-1) {
        psi1[j] =  psi0[j] - 0.5*B_B*(phi0[j-1] - 2.0*phi0[j] + phi0[j+1]); //+ 0.5*uj[j]*phi0[j];
        phi1[j] =  phi0[j] + 0.5*B_B*(psi0[j-1] - 2.0*psi0[j] + psi0[j+1]); //- 0.5*uj[j]*psi0[j];
    }
    psi1[NX-1] =  psi0[NX-1] - 0.5*B_B*(phi0[NX-2] - 2.0*phi0[NX-1]); //+ 0.5*uj[NX-1]*phi0[NX-1];
    phi1[NX-1] =  phi0[NX-1] + 0.5*B_B*(psi0[NX-2] - 2.0*psi0[NX-1]); //- 0.5*uj[NX-1]*psi0[NX-1];

    // Write arrays to file
    let filepath = format!("{}/{}-{}.dat",dirname,"wavepacket_snapshot",0);
    fs::File::create(&filepath).expect("Error creating file");
    let mut myfile = fs::OpenOptions::new()
        .write(true)
        .append(true)
        .open(filepath)
        .unwrap();

    let mut iline: String;

    for j in 0..NX {
        iline = format!("{} {} {} {}",x[j], uj[j], psi1[j], phi1[j]);
        writeln!(myfile, "{}", iline).expect("Error writing to file");
    }

    psi0 = psi1.clone();
    phi0 = phi1.clone();

    // Propagation in time

    for n in 1..NT {

        psi2[0] =  psi0[0] - B_B*(- 2.0*phi1[0] + phi1[1]) + uj[0]*phi1[0];
        phi2[0] =  phi0[0] + B_B*(- 2.0*psi1[0] + psi1[1]) - uj[0]*psi1[0];   
        for j in 1..(NX-1) {
            psi2[j] =  psi0[j] - B_B*(phi1[j-1] - 2.0*phi1[j] + phi1[j+1]) + uj[j]*phi1[j];
            phi2[j] =  phi0[j] + B_B*(psi1[j-1] - 2.0*psi1[j] + psi1[j+1]) - uj[j]*psi1[j];
        }
        psi2[NX-1] =  psi0[NX-1] - B_B*(phi1[NX-2] - 2.0*phi1[NX-1]) + uj[NX-1]*phi1[NX-1];
        phi2[NX-1] =  phi0[NX-1] + B_B*(psi1[NX-2] - 2.0*psi1[NX-1]) - uj[NX-1]*psi1[NX-1];
    
        // Write arrays to file
        if n as i32 % ST == 0 {
            let filepath = format!("{}/{}-{}.dat",dirname,"wavepacket_snapshot",n);
            fs::File::create(&filepath).expect("Error creating file");
            let mut myfile = fs::OpenOptions::new()
                .write(true)
                .append(true)
                .open(filepath)
                .unwrap();
        
            let mut iline: String;
        
            for j in 0..NX {
                iline = format!("{} {} {} {}",x[j], uj[j], psi2[j], phi2[j]);
                writeln!(myfile, "{}", iline).expect("Error writing to file");
            }
        }

        psi0 = psi1.clone();
        phi0 = phi1.clone();
        psi1 = psi2.clone();
        phi1 = phi2.clone();

    }

}
