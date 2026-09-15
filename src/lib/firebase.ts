import { initializeApp, getApps, getApp } from 'firebase/app';
import { getAuth } from 'firebase/auth';
import { getFirestore } from 'firebase/firestore';
import firebaseConfigJson from '../../firebase-applet-config.json';

// Allow environment variables to override the JSON config for easier deployment
const firebaseConfig = {
  apiKey: (typeof process !== 'undefined' && process.env?.VITE_FIREBASE_API_KEY) || firebaseConfigJson.apiKey,
  authDomain: (typeof process !== 'undefined' && process.env?.VITE_FIREBASE_AUTH_DOMAIN) || firebaseConfigJson.authDomain,
  projectId: (typeof process !== 'undefined' && process.env?.VITE_FIREBASE_PROJECT_ID) || firebaseConfigJson.projectId,
  storageBucket: (typeof process !== 'undefined' && process.env?.VITE_FIREBASE_STORAGE_BUCKET) || firebaseConfigJson.storageBucket,
  messagingSenderId: (typeof process !== 'undefined' && process.env?.VITE_FIREBASE_MESSAGING_SENDER_ID) || firebaseConfigJson.messagingSenderId,
  appId: (typeof process !== 'undefined' && process.env?.VITE_FIREBASE_APP_ID) || firebaseConfigJson.appId,
  firestoreDatabaseId: (typeof process !== 'undefined' && process.env?.VITE_FIREBASE_FIRESTORE_DATABASE_ID) || firebaseConfigJson.firestoreDatabaseId || '(default)'
};

const app = getApps().length === 0 ? initializeApp(firebaseConfig) : getApp();
export const db = getFirestore(app, firebaseConfig.firestoreDatabaseId);
export const auth = getAuth(app);
export default app;
